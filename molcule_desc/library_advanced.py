"""Higher-level reusable fragments composed from basic library structures."""

from .core import Structure, single, substitute
from .library import methane, phenyl


def _tert_butyl() -> Structure:
    """Build tert-butyl from one central methane and three methyl guests.

    Only the fourth port on the central carbon is exposed as the attachment
    point. All unused ports on the three guest methanes are capped with H.
    """
    return substitute(
        methane,
        {
            single(1): (methane, single(1)),
            single(2): (methane, single(1)),
            single(3): (methane, single(1)),
        },
        expose={single(1): (0, single(4))},
        name="tert-butyl",
    )


tert_butyl = _tert_butyl()


def _biphenyl() -> Structure:
    """Build biphenyl by connecting two phenyl fragments."""
    return substitute(
        phenyl,
        {single(1): (phenyl, single(1))},
        name="biphenyl",
    )


biphenyl = _biphenyl()


# BINOL-family chiral biaryl scaffolds.
# The two naphthalene units use the same Kekule topology as library.naphthalene.
# Atom 0/10 form the 1,1'-biaryl bond; atoms 1/11 carry the 2,2'-hydroxy groups.
from .core import Bond, BondType
from .library import _graph_structure


def _binol(name: str, substituents: dict[int, str] | None = None) -> Structure:
    """Build a BINOL scaffold, optionally substituting numbered ring carbons.

    Substituent keys are graph atom indices: 2/12 are the 3,3'-positions,
    while 8/18 are the 6,6'-positions in this topology.
    Supported substituents are Br, I, and methyl (C).
    """
    substituents = substituents or {}
    elements = ["C"] * 20 + ["O", "O"]
    edges = [
        (0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
        (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
        (4, 5, BondType.DOUBLE), (5, 0, BondType.SINGLE),
        (5, 6, BondType.SINGLE), (6, 7, BondType.DOUBLE),
        (7, 8, BondType.SINGLE), (8, 9, BondType.DOUBLE),
        (9, 4, BondType.SINGLE),
        (10, 11, BondType.DOUBLE), (11, 12, BondType.SINGLE),
        (12, 13, BondType.DOUBLE), (13, 14, BondType.SINGLE),
        (14, 15, BondType.DOUBLE), (15, 10, BondType.SINGLE),
        (15, 16, BondType.SINGLE), (16, 17, BondType.DOUBLE),
        (17, 18, BondType.SINGLE), (18, 19, BondType.DOUBLE),
        (19, 14, BondType.SINGLE),
        (0, 10, BondType.SINGLE),
        (1, 20, BondType.SINGLE), (11, 21, BondType.SINGLE),
    ]
    for atom_id, substituent in substituents.items():
        if atom_id not in range(20):
            raise ValueError(f"Invalid BINOL ring atom index: {atom_id}")
        if substituent not in {"Br", "I", "C"}:
            raise ValueError(f"Unsupported BINOL substituent: {substituent}")
        guest_id = len(elements)
        elements.append(substituent)
        edges.append((atom_id, guest_id, BondType.SINGLE))
    return _graph_structure(name, elements, edges)


binol = _binol("BINOL")
binol_3_3_dibromo = _binol("3,3'-dibromo-BINOL", {2: "Br", 12: "Br"})
binol_3_3_diiodo = _binol("3,3'-diiodo-BINOL", {2: "I", 12: "I"})
binol_6_6_dibromo = _binol("6,6'-dibromo-BINOL", {8: "Br", 18: "Br"})
binol_3_3_dimethyl = _binol("3,3'-dimethyl-BINOL", {2: "C", 12: "C"})
