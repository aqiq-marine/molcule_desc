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
    Supported substituents are Br, I, methyl (C), tert-butyl (tBu), and phenyl.
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
        if substituent not in {"Br", "I", "C", "tBu", "phenyl"}:
            raise ValueError(f"Unsupported BINOL substituent: {substituent}")
        guest_id = len(elements)
        if substituent in {"Br", "I", "C"}:
            elements.append(substituent)
            edges.append((atom_id, guest_id, BondType.SINGLE))
        elif substituent == "tBu":
            # Attach the tertiary carbon, then its three methyl groups.
            elements.extend(["C", "C", "C", "C"])
            edges.append((atom_id, guest_id, BondType.SINGLE))
            edges.extend(
                (guest_id, guest_id + offset, BondType.SINGLE)
                for offset in (1, 2, 3)
            )
        else:
            # A phenyl substituent with an alternating Kekule ring.
            elements.extend(["C"] * 6)
            edges.append((atom_id, guest_id, BondType.SINGLE))
            edges.extend([
                (guest_id, guest_id + 1, BondType.DOUBLE),
                (guest_id + 1, guest_id + 2, BondType.SINGLE),
                (guest_id + 2, guest_id + 3, BondType.DOUBLE),
                (guest_id + 3, guest_id + 4, BondType.SINGLE),
                (guest_id + 4, guest_id + 5, BondType.DOUBLE),
                (guest_id + 5, guest_id, BondType.SINGLE),
            ])
    return _graph_structure(name, elements, edges)


binol = _binol("BINOL")
binol_3_3_dibromo = _binol("3,3'-dibromo-BINOL", {2: "Br", 12: "Br"})
binol_3_3_diiodo = _binol("3,3'-diiodo-BINOL", {2: "I", 12: "I"})
binol_6_6_dibromo = _binol("6,6'-dibromo-BINOL", {8: "Br", 18: "Br"})
binol_3_3_dimethyl = _binol("3,3'-dimethyl-BINOL", {2: "C", 12: "C"})



def _binol_ring_edges() -> list[tuple[int, int, BondType]]:
    """Return the two fused naphthalene Kekule graphs joined at 1,1'."""
    return [
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
    ]


def _binol_diformyl() -> Structure:
    """Build 2,2'-diformyl-1,1'-binaphthyl (the OH groups are replaced)."""
    elements = ["C"] * 22 + ["O", "O"]
    edges = _binol_ring_edges()
    # Ring atoms 1 and 11 are the 2,2'-positions. Each now bears -CHO.
    edges.extend([
        (1, 20, BondType.SINGLE), (20, 22, BondType.DOUBLE),
        (11, 21, BondType.SINGLE), (21, 23, BondType.DOUBLE),
    ])
    return _graph_structure("2,2'-diformyl-1,1'-binaphthyl", elements, edges)


binol_2_2_diformyl = _binol_diformyl()


def _binol_bis_imine(amine: str) -> Structure:
    """Condense both formyl groups with an aminopyridine or aminopyrimidine.

    The imine N is connected to the exocyclic amino position (C2) of each
    heteroaromatic ring. Supported amines: 2-aminopyridine and
    2-aminopyrimidine.
    """
    if amine not in {"2-aminopyridine", "2-aminopyrimidine"}:
        raise ValueError(f"Unsupported imine partner: {amine}")
    elements = ["C"] * 22 + ["N", "N"]
    edges = _binol_ring_edges()
    # Formyl carbon atoms 20 and 21 form C=N imines; the carbonyl O is absent.
    edges.extend([
        (1, 20, BondType.SINGLE), (20, 22, BondType.DOUBLE),
        (11, 21, BondType.SINGLE), (21, 23, BondType.DOUBLE),
    ])
    for side, (formyl_c, imine_n) in enumerate(((20, 22), (21, 23))):
        start = len(elements)
        ring_elements = ["N", "C", "C", "C", "C", "C"]
        if amine == "2-aminopyrimidine":
            ring_elements = ["N", "C", "N", "C", "C", "C"]
        elements.extend(ring_elements)
        # The imine N is bonded to C2 of the heteroaromatic ring.
        edges.extend([
            (imine_n, start + 1, BondType.SINGLE),
            (start, start + 1, BondType.SINGLE),
            (start + 1, start + 2, BondType.DOUBLE),
            (start + 2, start + 3, BondType.SINGLE),
            (start + 3, start + 4, BondType.DOUBLE),
            (start + 4, start + 5, BondType.SINGLE),
            (start + 5, start, BondType.DOUBLE),
        ])
        if amine == "2-aminopyrimidine":
            # Replace the corresponding carbon in the ring with a second N.
            # This position has two single bonds in the Kekule form.
            edges[-5] = (start + 2, start + 3, BondType.SINGLE)
            edges[-4] = (start + 3, start + 4, BondType.DOUBLE)
            edges[-3] = (start + 4, start + 5, BondType.SINGLE)
            edges[-2] = (start + 5, start, BondType.DOUBLE)
    return _graph_structure(f"BINOL bis({amine} imine)", elements, edges)


binol_bis_2_pyridyl_imine = _binol_bis_imine("2-aminopyridine")
binol_bis_2_pyrimidyl_imine = _binol_bis_imine("2-aminopyrimidine")


# Friedel-Crafts-style substitution at the 6,6'-positions (graph atoms 8/18).
# These are structural examples; they do not encode a synthetic procedure.
binol_6_6_di_tert_butyl = _binol(
    "6,6'-di-tert-butyl-BINOL", {8: "tBu", 18: "tBu"}
)
binol_6_6_diphenyl = _binol(
    "6,6'-diphenyl-BINOL", {8: "phenyl", 18: "phenyl"}
)
