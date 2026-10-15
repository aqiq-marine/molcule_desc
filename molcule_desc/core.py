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



def substitute(
    host: Structure,
    substitutions: Mapping[Port, tuple[Structure, Port]],
    *,
    name: str | None = None,
) -> Structure:
    """Replace multiple exposed host ports with guest structures.

    substitutions maps each original host port to a (guest, guest_port) pair.
    Each guest is connected through its selected port. Unused ports on the
    host and guests remain exposed across intermediate substitutions.

    Host port keys always refer to the original host structure.
    """
    if not substitutions:
        return Structure(name or host.name, host.atoms.values(), host.bonds, host.ports)

    missing = set(substitutions) - set(host.ports)
    if missing:
        raise ValueError(f"Host has no exposed ports: {sorted(missing, key=repr)}")

    current = host
    host_port_map: dict[Port, Port] = {port: port for port in host.ports}
    guest_port_map: list[Port] = []

    for step, (original_port, (guest, guest_port)) in enumerate(
        substitutions.items(), start=1
    ):
        if original_port not in host_port_map:
            raise ValueError(f"Host port {original_port} has already been substituted")
        current_port = host_port_map[original_port]
        if guest_port not in guest.ports:
            raise ValueError(f"{guest.name} has no exposed port {guest_port}")

        next_index: dict[BondType, int] = {}
        expose: dict[Port, PortRef] = {}
        new_host_port_map: dict[Port, Port] = {}
        new_guest_port_map: list[Port] = []

        def publish(source: int, old_port: Port) -> Port:
            index = next_index.get(old_port.bond, 1)
            next_index[old_port.bond] = index + 1
            new_port = Port(old_port.bond, index)
            expose[new_port] = PortRef(source, old_port)
            return new_port

        for original, old_port in host_port_map.items():
            if original == original_port:
                continue
            new_host_port_map[original] = publish(0, old_port)

        for old_port in guest_port_map:
            new_guest_port_map.append(publish(0, old_port))

        for old_port in guest.ports:
            if old_port != guest_port:
                new_guest_port_map.append(publish(1, old_port))

        current = bond(
            current,
            current_port,
            guest,
            guest_port,
            expose=expose,
            name=name if step == len(substitutions) else f"{host.name} substitution {step}",
        )
        host_port_map = new_host_port_map
        guest_port_map = new_guest_port_map

    current.name = name or f"{host.name} substituted"
    return current


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
