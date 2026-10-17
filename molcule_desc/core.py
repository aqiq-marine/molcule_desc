"""Core graph and port operations for compositional molecular structures.

This prototype deliberately performs port and valence checks at runtime.
A Port is an exposed valence on an atom; connecting two ports consumes both
port declarations and creates a bond between their anchor atoms.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Mapping


class BondType(Enum):
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3


@dataclass(frozen=True)
class Port:
    bond: BondType
    index: int

    def __post_init__(self) -> None:
        if self.index < 1:
            raise ValueError("Port index must be a positive integer")


@dataclass(frozen=True)
class Atom:
    id: int
    element: str
    valence: int
    formal_charge: int = 0


@dataclass(frozen=True)
class Bond:
    start: int
    end: int
    order: BondType


@dataclass(frozen=True)
class PortRef:
    """Reference to a port on an input to bond(): 0=left, 1=right."""

    source: int
    port: Port

    def __post_init__(self) -> None:
        if self.source not in (0, 1):
            raise ValueError("PortRef.source must be 0 (left) or 1 (right)")


_VALENCE_LIMITS = {
    "H": 1,
    "B": 3,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "P": 5,
    "S": 6,
}


class Structure:
    """A molecular graph plus its currently exposed ports."""

    def __init__(
        self,
        name: str,
        atoms: Iterable[Atom],
        bonds: Iterable[Bond] = (),
        ports: Mapping[Port, int] | None = None,
    ) -> None:
        self.name = name
        atom_list = list(atoms)
        self.atoms = {atom.id: atom for atom in atom_list}
        if len(self.atoms) == 0:
            raise ValueError("A structure must contain at least one atom")
        if len(self.atoms) != len(atom_list):
            raise ValueError("Atom IDs must be unique")

        self.bonds = list(bonds)
        self.ports = dict(ports or {})

        for bond in self.bonds:
            if bond.start not in self.atoms or bond.end not in self.atoms:
                raise ValueError("Bond refers to an atom that does not exist")
            if bond.start == bond.end:
                raise ValueError("Self-bonds are not supported")

        for port, atom_id in self.ports.items():
            if atom_id not in self.atoms:
                raise ValueError(f"Port {port} refers to a missing atom")

        self.validate_valence()

    def validate_valence(self) -> None:
        """Check known elements against bond orders plus exposed-port orders."""
        valence = {atom_id: 0 for atom_id in self.atoms}
        for bond in self.bonds:
            order = bond.order.value
            valence[bond.start] += order
            valence[bond.end] += order
        for port, atom_id in self.ports.items():
            valence[atom_id] += port.bond.value

        for atom_id, used in valence.items():
            atom = self.atoms[atom_id]
            limit = atom.valence
            if limit is not None and used > limit:
                raise ValueError(
                    f"Valence exceeded for {atom.element} atom {atom_id}: "
                    f"{used} > {limit}"
                )

    @property
    def formula(self) -> dict[str, int]:
        """Return the explicit-atom formula. Call finalize() to cap ports."""
        return dict(Counter(atom.element for atom in self.atoms.values()))

    def formula_string(self) -> str:
        counts = self.formula
        order = ["C", "H"] + sorted(e for e in counts if e not in {"C", "H"})
        return "".join(
            element + (str(counts[element]) if counts[element] != 1 else "")
            for element in order
            if counts.get(element, 0)
        )

    def finalize(self) -> "Structure":
        """Resolve every exposed port to its terminal atom/group.

        A single port adds H, except for an H atom's own port (XS), which
        simply becomes that terminal hydrogen. A double port adds CH2; a
        triple port adds CH. The returned structure has no exposed ports.
        """
        atoms = dict(self.atoms)
        bonds = list(self.bonds)
        next_id = max(atoms, default=-1) + 1

        def add_atom(element: str) -> int:
            nonlocal next_id
            atom_id = next_id
            next_id += 1
            atoms[atom_id] = Atom(atom_id, element, _VALENCE_LIMITS[element])
            return atom_id

        for port, anchor in self.ports.items():
            if port.bond is BondType.SINGLE:
                if atoms[anchor].element == "H":
                    # XS: the terminal atom is already hydrogen.
                    continue
                hydrogen = add_atom("H")
                bonds.append(Bond(anchor, hydrogen, BondType.SINGLE))
            else:
                terminal_carbon = add_atom("C")
                bonds.append(Bond(anchor, terminal_carbon, port.bond))
                hydrogen_count = 2 if port.bond is BondType.DOUBLE else 1
                for _ in range(hydrogen_count):
                    hydrogen = add_atom("H")
                    bonds.append(
                        Bond(terminal_carbon, hydrogen, BondType.SINGLE)
                    )

        return Structure(
            f"{self.name} (finalized)",
            atoms.values(),
            bonds,
            {},
        )

    def __repr__(self) -> str:
        return (
            f"Structure(name={self.name!r}, atoms={len(self.atoms)}, "
            f"bonds={len(self.bonds)}, ports={list(self.ports)})"
        )


def atom(
    element: str,
    valence: int,
    *,
    formal_charge: int = 0,
) -> Atom:
    """Create an atom specification with an explicitly declared valence."""
    if valence < 0:
        raise ValueError("Valence must be non-negative")
    return Atom(0, element, valence, formal_charge)


def ported_atom(
    atom: Atom,
    ports: Iterable[Port],
    *,
    name: str | None = None,
) -> Structure:
    """Create a one-atom Structure from an Atom and explicit attachment ports."""
    port_list = list(ports)
    if len(port_list) != len(set(port_list)):
        raise ValueError("Duplicate port identifiers")
    return Structure(
        name or f"{atom.element} template",
        [atom],
        ports={port: atom.id for port in port_list},
    )


def single(index: int) -> Port:
    return Port(BondType.SINGLE, index)


def double(index: int) -> Port:
    return Port(BondType.DOUBLE, index)


def triple(index: int) -> Port:
    return Port(BondType.TRIPLE, index)


def bond(
    left: Structure,
    left_port: Port,
    right: Structure,
    right_port: Port,
    *,
    expose: Mapping[Port, PortRef],
    name: str | None = None,
) -> Structure:
    """Connect compatible ports on two structures.

    Each input structure is copied into the result. Use ring() to connect
    two ports within one existing structure without duplicating its graph.
    """
    return _connect(
        left, left_port, right, right_port,
        expose=expose, name=name, intramolecular=False,
    )


def _connect(
    left: Structure,
    left_port: Port,
    right: Structure,
    right_port: Port,
    *,
    expose: Mapping[Port, PortRef],
    name: str | None = None,
    intramolecular: bool,
) -> Structure:
    """Shared implementation for intermolecular composition and ring closure.

    expose maps each output-local Port to the input port it represents.
    PortRef.source is 0 for left and 1 for right.
    """
    if intramolecular and left is not right:
        raise ValueError("Intramolecular bonding requires the same Structure object")
    if left_port not in left.ports:
        raise ValueError(f"{left.name} has no exposed port {left_port}")
    if right_port not in right.ports:
        raise ValueError(f"{right.name} has no exposed port {right_port}")
    if left_port.bond is not right_port.bond:
        raise ValueError("Bond types do not match")
    if intramolecular and left_port == right_port:
        raise ValueError("A port cannot be connected to itself")

    same_structure = intramolecular
    atoms: dict[int, Atom] = {}
    bonds: list[Bond] = []
    maps: list[dict[int, int]] = [{}, {}]

    def copy_input(structure: Structure, side: int, offset: int) -> int:
        for atom_id, atom in structure.atoms.items():
            new_id = atom_id + offset
            maps[side][atom_id] = new_id
            atoms[new_id] = Atom(new_id, atom.element, atom.valence, atom.formal_charge)
        for edge in structure.bonds:
            bonds.append(
                Bond(
                    maps[side][edge.start],
                    maps[side][edge.end],
                    edge.order,
                )
            )
        return max(maps[side].values(), default=-1) + 1

    next_id = copy_input(left, 0, 0)
    if same_structure:
        maps[1] = maps[0]
    else:
        copy_input(right, 1, next_id)

    left_anchor = maps[0][left.ports[left_port]]
    right_anchor = maps[1][right.ports[right_port]]
    if left_anchor == right_anchor:
        raise ValueError("Cannot connect two ports anchored on the same atom")

    bonds.append(Bond(left_anchor, right_anchor, left_port.bond))

    input_ports = [left.ports, right.ports]
    consumed = [{left_port}, {right_port}]
    if same_structure:
        consumed[0].add(right_port)
    output_ports: dict[Port, int] = {}
    published_sources: set[tuple[int, Port]] = set()

    for output_port, ref in expose.items():
        source_index = 0 if same_structure else ref.source
        source_ports = input_ports[source_index]
        if ref.port not in source_ports:
            raise ValueError(f"Input has no exposed port {ref.port}")
        if ref.port in consumed[source_index]:
            raise ValueError(f"Cannot expose consumed port {ref.port}")
        if output_port.bond is not ref.port.bond:
            raise ValueError("Output port bond type does not match source port")
        source_key = (source_index, ref.port)
        if source_key in published_sources:
            raise ValueError(f"Input port {ref.port} is published more than once")
        published_sources.add(source_key)
        source_atom = source_ports[ref.port]
        output_ports[output_port] = maps[source_index][source_atom]

    # Any unconsumed input port not explicitly exported is finalized locally.
    # This preserves the rule that only exported ports remain available.
    next_id = max(atoms, default=-1) + 1

    def cap_unexported(port: Port, anchor: int) -> None:
        nonlocal next_id
        if port.bond is BondType.SINGLE:
            if atoms[anchor].element == "H":
                return
            terminal_id = next_id
            next_id += 1
            atoms[terminal_id] = Atom(terminal_id, "H", 1)
            bonds.append(Bond(anchor, terminal_id, BondType.SINGLE))
            return

        terminal_id = next_id
        next_id += 1
        atoms[terminal_id] = Atom(terminal_id, "C", _VALENCE_LIMITS["C"])
        bonds.append(Bond(anchor, terminal_id, port.bond))
        hydrogen_count = 2 if port.bond is BondType.DOUBLE else 1
        for _ in range(hydrogen_count):
            hydrogen_id = next_id
            next_id += 1
            atoms[hydrogen_id] = Atom(hydrogen_id, "H", 1)
            bonds.append(Bond(terminal_id, hydrogen_id, BondType.SINGLE))

    source_count = 1 if same_structure else 2
    for source_index in range(source_count):
        for port, atom_id in input_ports[source_index].items():
            if port in consumed[source_index]:
                continue
            if (source_index, port) in published_sources:
                continue
            cap_unexported(port, maps[source_index][atom_id])

    return Structure(
        name or f"bond({left.name}, {right.name})",
        atoms.values(),
        bonds,
        output_ports,
    )




def fuse(
    cyclic: Structure,
    chain: Structure,
    *,
    cyclic_ports: tuple[Port, Port],
    chain_ports: tuple[Port, Port],
    name: str | None = None,
) -> Structure:
    """Fuse a chain to a cyclic structure by consuming two port pairs.

    The first port in each tuple is connected first; the second pair closes
    the new ring. Every other exposed port from both input structures remains
    exposed on the returned Structure, so the fused structure stays editable.
    Port labels in the result are reassigned uniquely by bond type.
    """
    if len(set(cyclic_ports)) != 2:
        raise ValueError("cyclic_ports must contain two distinct ports")
    if len(set(chain_ports)) != 2:
        raise ValueError("chain_ports must contain two distinct ports")
    for port in cyclic_ports:
        if port not in cyclic.ports:
            raise ValueError(f"{cyclic.name} has no exposed port {port}")
    for port in chain_ports:
        if port not in chain.ports:
            raise ValueError(f"{chain.name} has no exposed port {port}")

    # bond() consumes the first pair. Export every other port, including the
    # second pair needed for ring closure, while assigning unique output labels.
    expose: dict[Port, PortRef] = {}
    source_to_output: dict[tuple[int, Port], Port] = {}
    next_index: dict[BondType, int] = {}
    for source, structure, consumed in (
        (0, cyclic, cyclic_ports[0]),
        (1, chain, chain_ports[0]),
    ):
        for port in structure.ports:
            if port == consumed:
                continue
            index = next_index.get(port.bond, 1)
            next_index[port.bond] = index + 1
            output_port = Port(port.bond, index)
            expose[output_port] = PortRef(source, port)
            source_to_output[(source, port)] = output_port

    attached = bond(
        cyclic,
        cyclic_ports[0],
        chain,
        chain_ports[0],
        expose=expose,
        name=f"{name or 'fused'} (open)",
    )
    closure_first = source_to_output[(0, cyclic_ports[1])]
    closure_second = source_to_output[(1, chain_ports[1])]
    return ring(
        attached,
        closure_first,
        closure_second,
        expose={port: port for port in attached.ports if port not in {
            closure_first, closure_second
        }},
        name=name or f"fuse({cyclic.name}, {chain.name})",
    )

def substitute(
    host: Structure,
    substitutions: Mapping[Port, tuple[Structure, Port]],
    *,
    expose: Mapping[Port, tuple[int, Port]] | None = None,
    name: str | None = None,
) -> Structure:
    """Replace multiple exposed host ports with guest structures.

    substitutions maps each original host port to a (guest, guest_port) pair.
    expose optionally maps each output port to (source, port), where source 0
    is the original host and sources 1..N are guests in substitutions order.

    If expose is omitted, all unused host and guest ports are preserved.
    If supplied, only listed ports are exposed; other unused ports are capped.
    """
    entries = list(substitutions.items())
    missing = {host_port for host_port, _ in entries} - set(host.ports)
    if missing:
        raise ValueError(f"Host has no exposed ports: {sorted(missing, key=repr)}")

    structures = [host] + [guest for _, (guest, _) in entries]
    consumed: list[set[Port]] = [set() for _ in structures]
    for host_port, (guest, guest_port) in entries:
        if guest_port not in guest.ports:
            raise ValueError(f"{guest.name} has no exposed port {guest_port}")
        if host_port.bond is not guest_port.bond:
            raise ValueError("Bond types do not match")
        consumed[0].add(host_port)
    for index, (_, (_, guest_port)) in enumerate(entries, start=1):
        consumed[index].add(guest_port)

    atoms: dict[int, Atom] = {}
    bonds: list[Bond] = []
    atom_maps: list[dict[int, int]] = []
    next_atom_id = 0
    for structure in structures:
        mapping: dict[int, int] = {}
        for atom_id, atom_spec in structure.atoms.items():
            mapping[atom_id] = next_atom_id
            atoms[next_atom_id] = Atom(
                next_atom_id, atom_spec.element, atom_spec.valence,
                atom_spec.formal_charge,
            )
            next_atom_id += 1
        atom_maps.append(mapping)
        for edge in structure.bonds:
            bonds.append(Bond(mapping[edge.start], mapping[edge.end], edge.order))

    for index, (host_port, (guest, guest_port)) in enumerate(entries, start=1):
        host_anchor = atom_maps[0][host.ports[host_port]]
        guest_anchor = atom_maps[index][guest.ports[guest_port]]
        if host_anchor == guest_anchor:
            raise ValueError("Cannot connect ports anchored on the same atom")
        bonds.append(Bond(host_anchor, guest_anchor, host_port.bond))

    available: dict[tuple[int, Port], int] = {}
    for source, structure in enumerate(structures):
        for port, atom_id in structure.ports.items():
            if port not in consumed[source]:
                available[(source, port)] = atom_maps[source][atom_id]

    if expose is None:
        output_ports: dict[Port, int] = {}
        next_index: dict[BondType, int] = {}
        for (source, port), anchor in available.items():
            index = next_index.get(port.bond, 1)
            next_index[port.bond] = index + 1
            output_ports[Port(port.bond, index)] = anchor
        exported_keys = set(available)
    else:
        output_ports = {}
        exported_keys: set[tuple[int, Port]] = set()
        for output_port, (source, source_port) in expose.items():
            if source < 0 or source >= len(structures):
                raise ValueError(
                    f"Invalid expose source {source}; expected 0..{len(structures) - 1}"
                )
            key = (source, source_port)
            if key not in available:
                if source_port in consumed[source]:
                    raise ValueError(f"Cannot expose consumed port {source_port}")
                raise ValueError(f"Input has no exposed port {source_port}")
            if output_port.bond is not source_port.bond:
                raise ValueError("Output port bond type does not match source port")
            if key in exported_keys:
                raise ValueError(f"Input port {source_port} is published more than once")
            exported_keys.add(key)
            output_ports[output_port] = available[key]

    next_atom_id = max(atoms, default=-1) + 1
    for key, anchor in available.items():
        if key in exported_keys:
            continue
        port = key[1]
        if port.bond is BondType.SINGLE:
            if atoms[anchor].element == "H":
                continue
            atoms[next_atom_id] = Atom(next_atom_id, "H", 1)
            bonds.append(Bond(anchor, next_atom_id, BondType.SINGLE))
            next_atom_id += 1
        else:
            terminal_id = next_atom_id
            atoms[terminal_id] = Atom(terminal_id, "C", _VALENCE_LIMITS["C"])
            bonds.append(Bond(anchor, terminal_id, port.bond))
            next_atom_id += 1
            hydrogen_count = 2 if port.bond is BondType.DOUBLE else 1
            for _ in range(hydrogen_count):
                atoms[next_atom_id] = Atom(next_atom_id, "H", 1)
                bonds.append(Bond(terminal_id, next_atom_id, BondType.SINGLE))
                next_atom_id += 1

    return Structure(name or f"{host.name} substituted", atoms.values(), bonds, output_ports)


def ring(
    structure: Structure,
    first_port: Port,
    second_port: Port,
    *,
    expose: Mapping[Port, Port],
    name: str | None = None,
) -> Structure:
    """Close two ports on one structure without duplicating its graph.

    The expose mapping assigns result-local ports directly to existing ports
    on the input structure, e.g. {single(1): single(3)}.
    """
    return _connect(
        structure,
        first_port,
        structure,
        second_port,
        expose={output_port: PortRef(0, source_port)
                for output_port, source_port in expose.items()},
        name=name or f"ring({structure.name})",
        intramolecular=True,
    )
