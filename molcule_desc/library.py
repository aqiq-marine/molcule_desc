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
    substitute,
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
    methyl = ported_atom(
        atom("C", 4),
        [single(1), single(2), single(3), single(4)],
        name="methyl C",
    )
    return substitute(
        _carbonyl(),
        {
            single(1): (methyl, single(1)),
            single(2): (methyl, single(1)),
        },
        expose={
            single(1): (1, single(2)),
            single(2): (1, single(3)),
            single(3): (1, single(4)),
            single(4): (2, single(2)),
            single(5): (2, single(3)),
            single(6): (2, single(4)),
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


# Allowed valences used when choosing a default valence for library atoms.
# These are valence/bonding capacities, not oxidation states. When a structure
# needs a particular state (e.g. sulfone sulfur), pass it explicitly via valences.
_ELEMENT_VALENCES: dict[str, tuple[int, ...]] = {
    "H": (1,),
    "B": (3,),
    "C": (4,),
    "N": (3, 5),
    "O": (2,),
    "F": (1,),
    "Cl": (1, 3, 5, 7),
    "Br": (1, 3, 5, 7),
    "I": (1, 3, 5, 7),
    "P": (3, 5),
    "S": (2, 4, 6),
}


def _graph_structure(
    name: str,
    elements: list[str],
    edges: list[tuple[int, int, BondType]],
    *,
    valences: list[int] | None = None,
    formal_charges: list[int] | None = None,
) -> Structure:
    """Build a graph and expose remaining valences as single-bond ports.

    For each atom, choose the smallest allowed valence that can accommodate
    its existing bond-order sum. Explicit valences override this inference
    when a particular bonding state is required by the fragment.
    """
    charges = formal_charges or [0] * len(elements)
    used = [0] * len(elements)
    bonds = [Bond(start, end, order) for start, end, order in edges]
    for edge in bonds:
        used[edge.start] += edge.order.value
        used[edge.end] += edge.order.value

    if valences is not None:
        if len(valences) != len(elements):
            raise ValueError("valences must have one entry per atom")
        atom_valences = list(valences)
    else:
        atom_valences = []
        for atom_id, element in enumerate(elements):
            allowed = _ELEMENT_VALENCES.get(element)
            if allowed is None:
                raise ValueError(f"No allowed valences configured for element {element!r}")
            candidates = [value for value in allowed if value >= used[atom_id]]
            if not candidates:
                raise ValueError(
                    f"Invalid library structure {name!r}: atom {atom_id} ({element}) "
                    f"uses bond-order sum {used[atom_id]}, but allowed valences are {allowed}"
                )
            atom_valences.append(min(candidates))

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

# A four-carbon chain with one exposed port at each terminal carbon.
butane = bond(
    propane,
    single(1),
    methane,
    single(1),
    expose={
        single(1): PortRef(0, single(6)),
        single(2): PortRef(1, single(2)),
    },
    name="butane",
)

def _linear_alkane(carbon_count: int) -> Structure:
    """Build a straight-chain alkane with every C-H valence exposed as a port."""
    if not 1 <= carbon_count <= 10:
        raise ValueError("carbon_count must be between 1 and 10")
    return _graph_structure(
        f"alkane_{carbon_count}",
        ["C"] * carbon_count,
        [(i, i + 1, BondType.SINGLE) for i in range(carbon_count - 1)],
    )


pentane = _linear_alkane(5)
hexane = _linear_alkane(6)
heptane = _linear_alkane(7)
octane = _linear_alkane(8)
nonane = _linear_alkane(9)
decane = _linear_alkane(10)


def _chain_preset(name: str, carbon_count: int) -> Structure:
    """Create a saturated carbon chain with every remaining valence exposed."""
    return _graph_structure(
        name,
        ["C"] * carbon_count,
        [(i, i + 1, BondType.SINGLE) for i in range(carbon_count - 1)],
    )


# Flexible chain fragments intended for ring fusion. All available valences,
# including those on internal carbons, remain exposed for later substitution.
trimethylene = _chain_preset("trimethylene", 3)
tetramethylene = _chain_preset("tetramethylene", 4)
pentamethylene = _chain_preset("pentamethylene", 5)
hexamethylene = _chain_preset("hexamethylene", 6)
heptamethylene = _chain_preset("heptamethylene", 7)
octamethylene = _chain_preset("octamethylene", 8)
nonamethylene = _chain_preset("nonamethylene", 9)
decamethylene = _chain_preset("decamethylene", 10)
acetylene = _graph_structure(
    "acetylene", ["C", "C"], [(0, 1, BondType.TRIPLE)]
)
# Cycloalkanes from three- to ten-membered rings.
# The ring closure consumes one valence on each terminal carbon; remaining
# valences are exposed as ports and become hydrogens on finalize().
def _cycloalkane(size: int) -> Structure:
    if not 3 <= size <= 10:
        raise ValueError("Cycloalkane ring size must be between 3 and 10")
    names = {
        3: "cyclopropane",
        4: "cyclobutane",
        5: "cyclopentane",
        6: "cyclohexane",
        7: "cycloheptane",
        8: "cyclooctane",
        9: "cyclononane",
        10: "cyclodecane",
    }
    return _graph_structure(
        names[size],
        ["C"] * size,
        _cycle_edges(list(range(size)), [BondType.SINGLE] * size),
    )


cyclopropane = _cycloalkane(3)
cyclobutane = _cycloalkane(4)
cyclopentane = _cycloalkane(5)
cyclohexane = _cycloalkane(6)
cycloheptane = _cycloalkane(7)
cyclooctane = _cycloalkane(8)
cyclononane = _cycloalkane(9)
cyclodecane = _cycloalkane(10)

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
pyrimidine = _graph_structure(
    "pyrimidine", ["N", "C", "N", "C", "C", "C"],
    _cycle_edges(
        [0, 1, 2, 3, 4, 5],
        [BondType.DOUBLE, BondType.SINGLE, BondType.DOUBLE,
         BondType.SINGLE, BondType.DOUBLE, BondType.SINGLE],
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
    valences=[6, 2, 2],
)
phosphine = _graph_structure("phosphine", ["P"], [], valences=[3])
phosphine_oxide = _graph_structure(
    "phosphine_oxide", ["P", "O"], [(0, 1, BondType.DOUBLE)],
    valences=[5, 2],
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


diphenyl_ether = _biaryl("diphenyl_ether", ether_link=True)


# Fused-ring aromatic compounds represented as Kekule graphs.
naphthol = _graph_structure(
    "naphthol", ["C"] * 10 + ["O"],
    [
        (0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
        (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
        (4, 5, BondType.DOUBLE), (5, 0, BondType.SINGLE),
        (5, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
        (7, 8, BondType.SINGLE), (8, 9, BondType.DOUBLE),
        (9, 4, BondType.SINGLE), (2, 10, BondType.SINGLE),
    ],
    valences=[4] * 10 + [2],
)

indole = _graph_structure(
    "indole", ["C"] * 8 + ["N"],
    [
        (0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
        (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
        (4, 5, BondType.DOUBLE), (5, 0, BondType.SINGLE),
        (4, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
        (7, 8, BondType.SINGLE), (8, 5, BondType.SINGLE),
    ],
    valences=[4] * 8 + [3],
)

benzofuran = _graph_structure(
    "benzofuran", ["C"] * 8 + ["O"],
    [
        (0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
        (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
        (4, 5, BondType.DOUBLE), (5, 0, BondType.SINGLE),
        (4, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
        (7, 8, BondType.SINGLE), (8, 5, BondType.SINGLE),
    ],
    valences=[4] * 8 + [2],
)

benzothiophene = _graph_structure(
    "benzothiophene", ["C"] * 8 + ["S"],
    [
        (0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
        (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
        (4, 5, BondType.DOUBLE), (5, 0, BondType.SINGLE),
        (4, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
        (7, 8, BondType.SINGLE), (8, 5, BondType.SINGLE),
    ],
    valences=[4] * 8 + [2],
)

quinoline = _graph_structure(
    "quinoline", ["N"] + ["C"] * 9,
    [
        (0, 1, BondType.SINGLE), (1, 2, BondType.DOUBLE),
        (2, 3, BondType.SINGLE), (3, 4, BondType.DOUBLE),
        (4, 5, BondType.SINGLE), (5, 0, BondType.DOUBLE),
        (5, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
        (7, 8, BondType.SINGLE), (8, 9, BondType.DOUBLE),
        (9, 4, BondType.SINGLE),
    ],
)

isoquinoline = _graph_structure(
    "isoquinoline", ["N"] + ["C"] * 9,
    [
        (0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
        (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
        (4, 5, BondType.DOUBLE), (5, 0, BondType.SINGLE),
        (5, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
        (7, 8, BondType.SINGLE), (8, 9, BondType.DOUBLE),
        (9, 4, BondType.SINGLE),
    ],
)

anthracene = _graph_structure(
    "anthracene", ["C"] * 14,
    [
        (0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
        (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
        (4, 5, BondType.DOUBLE), (5, 0, BondType.SINGLE),
        (5, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
        (7, 8, BondType.SINGLE), (8, 9, BondType.DOUBLE),
        (9, 4, BondType.SINGLE), (8, 10, BondType.SINGLE),
        (10, 11, BondType.DOUBLE), (11, 12, BondType.SINGLE),
        (12, 13, BondType.DOUBLE), (13, 9, BondType.SINGLE),
    ],
)
