"""Reusable elementary structures and example compounds."""

from .core import (
    Atom,
    Bond,
    BondType,
    Port,
    PortRef,
    Structure,
    atom,
    ported_atom,
    bond,
    double,
    single,
    triple,
    ring,
)

# These are atom-topology templates, not virtual atoms. Their names describe
# the bond ports available on each real atom.
hydrogen_single = ported_atom(atom("H", 1), [single(1)], name="hydrogen_single")
carbon_double_two_single = ported_atom(
    atom("C", 4), [double(1), single(1), single(2)],
    name="carbon_double_two_single",
)
carbon_triple_single = ported_atom(
    atom("C", 4), [triple(1), single(1)], name="carbon_triple_single"
)
oxygen_double = ported_atom(atom("O", 2), [double(1)], name="oxygen_double")


def _ethylene() -> Structure:
    return bond(
        carbon_double_two_single, double(1), carbon_double_two_single, double(1),
        expose={
            single(1): PortRef(0, single(1)),
            single(2): PortRef(0, single(2)),
            single(3): PortRef(1, single(1)),
            single(4): PortRef(1, single(2)),
        },
        name="ethylene",
    )


def _methanimine() -> Structure:
    carbon = ported_atom(atom("C", 4), [double(1), single(1), single(2)], name="imine C")
    nitrogen = ported_atom(atom("N", 3), [double(1), single(1)], name="imine N")
    return bond(
        carbon, double(1), nitrogen, double(1),
        expose={
            single(1): PortRef(0, single(1)),
            single(2): PortRef(0, single(2)),
            single(3): PortRef(1, single(1)),
        },
        name="methanimine",
    )


def _carbonyl() -> Structure:
    # carbon_double_two_single's double port is consumed to form C=O; its two single ports are
    # still exposed so the carbonyl carbon can be connected to substituents.
    return bond(
        carbon_double_two_single, double(1), oxygen_double, double(1),
        expose={
            single(1): PortRef(0, single(1)),
            single(2): PortRef(0, single(2)),
        },
        name="carbonyl",
    )


def _acetone() -> Structure:
    carbonyl = _carbonyl()
    methyl = ported_atom(
        atom("C", 4), [single(1), single(2), single(3), single(4)], name="methyl C"
    )
    left_attached = bond(
        carbonyl, single(1), methyl, single(1),
        expose={
            single(1): PortRef(0, single(2)),
            single(2): PortRef(1, single(2)),
            single(3): PortRef(1, single(3)),
            single(4): PortRef(1, single(4)),
        },
        name="acetone (first methyl attached)",
    )
    return bond(
        left_attached, single(1), methyl, single(1),
        expose={
            single(1): PortRef(0, single(2)),
            single(2): PortRef(0, single(3)),
            single(3): PortRef(0, single(4)),
            single(4): PortRef(1, single(2)),
            single(5): PortRef(1, single(3)),
            single(6): PortRef(1, single(4)),
        },
        name="acetone",
    )


ethylene = _ethylene()
methanimine = _methanimine()
carbonyl = _carbonyl()
formaldehyde = carbonyl
acetone = _acetone()


def _benzene() -> Structure:
    # Kekule representation with one exposed single-bond port for each C-H bond.
    # Ports 1 and 2 close the ring; ports 3-8 represent the six hydrogens.
    ports = {single(1): 0, single(2): 5}
    ports.update({single(i + 3): i for i in range(6)})
    chain = Structure(
        "benzene chain",
        [Atom(i, "C", 4) for i in range(6)],
        [
            Bond(0, 1, BondType.DOUBLE),
            Bond(1, 2, BondType.SINGLE),
            Bond(2, 3, BondType.DOUBLE),
            Bond(3, 4, BondType.SINGLE),
            Bond(4, 5, BondType.DOUBLE),
        ],
        ports=ports,
    )
    return ring(
        chain,
        single(1),
        single(2),
        expose={single(i + 1): single(i + 3) for i in range(6)},
        name="benzene",
    )


benzene = _benzene()


# Generic graph builder for reusable library structures. Every unused valence is
# represented by a single-bond port, so it can either be connected to a fragment
# or completed with hydrogen when the structure is finalized.
_ELEMENT_VALENCES = {
    "H": 1, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1,
    "Cl": 1, "Br": 1, "I": 1, "P": 5, "S": 6,
}


def _graph_structure(
    name: str,
    elements: list[str],
    edges: list[tuple[int, int, BondType]],
    *,
    valences: list[int] | None = None,
    formal_charges: list[int] | None = None,
) -> Structure:
    """Build a graph and expose every remaining valence as a single port."""
    atom_valences = valences or [_ELEMENT_VALENCES[element] for element in elements]
    charges = formal_charges or [0] * len(elements)
    used = [0] * len(elements)
    bonds = [Bond(start, end, order) for start, end, order in edges]
    for edge in bonds:
        used[edge.start] += edge.order.value
        used[edge.end] += edge.order.value
    atoms = [
        Atom(i, element, atom_valences[i], charges[i])
        for i, element in enumerate(elements)
    ]
    ports: dict[Port, int] = {}
    port_index = 1
    for atom_id, valence in enumerate(atom_valences):
        remaining = valence - used[atom_id]
        if remaining < 0:
            raise ValueError(
                f"Invalid library structure {name!r}: atom {atom_id} "
                f"uses valence {used[atom_id]} > {valence}"
            )
        for _ in range(remaining):
            ports[single(port_index)] = atom_id
            port_index += 1
    return Structure(name, atoms, bonds, ports)


def _cycle_edges(
    atom_ids: list[int],
    orders: list[BondType],
) -> list[tuple[int, int, BondType]]:
    return [
        (atom_ids[i], atom_ids[(i + 1) % len(atom_ids)], orders[i])
        for i in range(len(atom_ids))
    ]


# Saturated and unsaturated hydrocarbon skeletons.
methane = _graph_structure("methane", ["C"], [])
ethane = _graph_structure("ethane", ["C", "C"], [(0, 1, BondType.SINGLE)])
propane = _graph_structure(
    "propane", ["C", "C", "C"],
    [(0, 1, BondType.SINGLE), (1, 2, BondType.SINGLE)],
)
acetylene = _graph_structure(
    "acetylene", ["C", "C"], [(0, 1, BondType.TRIPLE)]
)
cyclohexane = _graph_structure(
    "cyclohexane", ["C"] * 6,
    _cycle_edges([0, 1, 2, 3, 4, 5], [BondType.SINGLE] * 6),
)

# These ring structures expose all substitutable positions. Consuming one port
# gives a phenyl substituent; consuming two gives a phenylene linker.
phenyl = _graph_structure(
    "phenyl", ["C"] * 6,
    _cycle_edges(
        [0, 1, 2, 3, 4, 5],
        [BondType.DOUBLE, BondType.SINGLE, BondType.DOUBLE,
         BondType.SINGLE, BondType.DOUBLE, BondType.SINGLE],
    ),
)
phenylene = _graph_structure(
    "phenylene", ["C"] * 6,
    _cycle_edges(
        [0, 1, 2, 3, 4, 5],
        [BondType.DOUBLE, BondType.SINGLE, BondType.DOUBLE,
         BondType.SINGLE, BondType.DOUBLE, BondType.SINGLE],
    ),
)

# Oxygen-containing groups and common carbonyl derivatives.
hydroxyl = _graph_structure("hydroxyl", ["O"], [])
ether = _graph_structure(
    "ether", ["C", "O", "C"],
    [(0, 1, BondType.SINGLE), (1, 2, BondType.SINGLE)],
)
carboxyl = _graph_structure(
    "carboxyl", ["C", "O", "O"],
    [(0, 1, BondType.DOUBLE), (0, 2, BondType.SINGLE)],
)
ester = _graph_structure(
    "ester", ["C", "O", "O", "C"],
    [(0, 1, BondType.DOUBLE), (0, 2, BondType.SINGLE),
     (2, 3, BondType.SINGLE)],
)

# Nitrogen-containing groups.
amino = _graph_structure("amino", ["N"], [])
amine = _graph_structure("amine", ["C", "N"], [(0, 1, BondType.SINGLE)])
amide = _graph_structure(
    "amide", ["C", "O", "N"],
    [(0, 1, BondType.DOUBLE), (0, 2, BondType.SINGLE)],
)
nitrile = _graph_structure(
    "nitrile", ["C", "N"], [(0, 1, BondType.TRIPLE)]
)
imine = _graph_structure("imine", ["C", "N"], [(0, 1, BondType.DOUBLE)])

# Aromatic rings are represented by Kekule graphs; RDKit perceives aromaticity
# when converting a completed structure to SMILES.
pyridine = _graph_structure(
    "pyridine", ["N", "C", "C", "C", "C", "C"],
    _cycle_edges(
        [0, 1, 2, 3, 4, 5],
        [BondType.SINGLE, BondType.DOUBLE, BondType.SINGLE,
         BondType.DOUBLE, BondType.SINGLE, BondType.DOUBLE],
    ),
)
pyrrole = _graph_structure(
    "pyrrole", ["N", "C", "C", "C", "C"],
    _cycle_edges(
        [0, 1, 2, 3, 4],
        [BondType.SINGLE, BondType.DOUBLE, BondType.SINGLE,
         BondType.DOUBLE, BondType.SINGLE],
    ),
)
thiophene = _graph_structure(
    "thiophene", ["S", "C", "C", "C", "C"],
    _cycle_edges(
        [0, 1, 2, 3, 4],
        [BondType.SINGLE, BondType.DOUBLE, BondType.SINGLE,
         BondType.DOUBLE, BondType.SINGLE],
    ),
    valences=[2, 4, 4, 4, 4],
)
imidazole = _graph_structure(
    "imidazole", ["N", "C", "N", "C", "C"],
    [(0, 1, BondType.SINGLE), (1, 2, BondType.DOUBLE),
     (2, 3, BondType.SINGLE), (3, 4, BondType.DOUBLE),
     (4, 0, BondType.SINGLE)],
)
naphthalene = _graph_structure(
    "naphthalene", ["C"] * 10,
    [(0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
     (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
     (4, 5, BondType.DOUBLE), (5, 0, BondType.SINGLE),
     (5, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
     (7, 8, BondType.SINGLE), (8, 9, BondType.DOUBLE),
     (9, 4, BondType.SINGLE)],
)

# Halogen substituents.
fluoro = _graph_structure("fluoro", ["F"], [])
chloro = _graph_structure("chloro", ["Cl"], [])
bromo = _graph_structure("bromo", ["Br"], [])
iodo = _graph_structure("iodo", ["I"], [])

# Sulfur- and phosphorus-containing groups.
thiol = _graph_structure("thiol", ["S"], [], valences=[2])
thioether = _graph_structure(
    "thioether", ["C", "S", "C"],
    [(0, 1, BondType.SINGLE), (1, 2, BondType.SINGLE)],
    valences=[4, 2, 4],
)
sulfone = _graph_structure(
    "sulfone", ["S", "O", "O"],
    [(0, 1, BondType.DOUBLE), (0, 2, BondType.DOUBLE)],
)
phosphine = _graph_structure("phosphine", ["P"], [], valences=[3])
phosphine_oxide = _graph_structure(
    "phosphine_oxide", ["P", "O"], [(0, 1, BondType.DOUBLE)]
)
phosphate_ester = _graph_structure(
    "phosphate_ester", ["P", "O", "O", "O", "O"],
    [(0, 1, BondType.DOUBLE), (0, 2, BondType.SINGLE),
     (0, 3, BondType.SINGLE), (0, 4, BondType.SINGLE)],
    valences=[5, 2, 2, 2, 2],
)

# Linkers and biaryl examples.
methylene = _graph_structure("methylene", ["C"], [])
ethanol = _graph_structure(
    "ethanol", ["C", "C", "O"],
    [(0, 1, BondType.SINGLE), (1, 2, BondType.SINGLE)],
)
ethylene_linker = _graph_structure(
    "ethylene_linker", ["C", "C"], [(0, 1, BondType.SINGLE)]
)


def _biaryl(name: str, *, ether_link: bool) -> Structure:
    elements = ["C"] * 12
    edges = (
        _cycle_edges(
            [0, 1, 2, 3, 4, 5],
            [BondType.DOUBLE, BondType.SINGLE, BondType.DOUBLE,
             BondType.SINGLE, BondType.DOUBLE, BondType.SINGLE],
        )
        + _cycle_edges(
            [6, 7, 8, 9, 10, 11],
            [BondType.DOUBLE, BondType.SINGLE, BondType.DOUBLE,
             BondType.SINGLE, BondType.DOUBLE, BondType.SINGLE],
        )
    )
    if ether_link:
        elements.append("O")
        edges.extend([(0, 12, BondType.SINGLE), (6, 12, BondType.SINGLE)])
    else:
        edges.append((0, 6, BondType.SINGLE))
    return _graph_structure(name, elements, edges)


biphenyl = _biaryl("biphenyl", ether_link=False)
diphenyl_ether = _biaryl("diphenyl_ether", ether_link=True)
