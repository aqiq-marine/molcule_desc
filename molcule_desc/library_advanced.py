"""Higher-level reusable fragments composed from basic library structures."""

from .core import Structure, single, substitute
from .library import (
    bromo, carbonyl, imine, iodo, methane, naphthalene, naphthol, phenyl, pyridine,
)


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


def _port_at_atom(structure: Structure, atom_id: int) -> Port:
    """Return an exposed single-bond port anchored to a specified atom."""
    for port, anchor in structure.ports.items():
        if anchor == atom_id:
            return port
    raise ValueError(f"{structure.name} has no exposed port on atom {atom_id}")


def _compose_binol_base() -> Structure:
    """Compose BINOL from two naphthol fragments joined at their 1-positions.

    In the naphthol fragment, atom 1 is the biaryl attachment position and
    atom 2 bears OH. The two fragments are joined through atom 1, leaving
    their remaining ports available for further substitution.
    """
    join_port = _port_at_atom(naphthol, 1)
    return substitute(
        naphthol,
        {join_port: (naphthol, join_port)},
        name="BINOL",
    )


def _replace_with_fragment(
    host: Structure,
    anchor_atom: int,
    guest: Structure,
    *,
    name: str,
) -> Structure:
    """Attach a guest fragment at a host atom and cap all guest-side spare ports."""
    host_port = _port_at_atom(host, anchor_atom)
    expose = {
        single(index): (0, port)
        for index, port in enumerate(
            (p for p in host.ports if p != host_port), start=1
        )
    }
    return substitute(
        host,
        {host_port: (guest, _port_at_atom(guest, next(iter(guest.ports.values()))))},
        expose=expose,
        name=name,
    )


def _binol(name: str, substituents: dict[int, str] | None = None) -> Structure:
    """Build BINOL and derivatives compositionally from reusable fragments.

    Substituent keys are atom indices in the composed two-naphthol scaffold.
    Supported substituents: Br, I, methyl (C), tert-butyl (tBu), and phenyl.
    """
    result = _compose_binol_base()
    guests = {
        "Br": bromo,
        "I": iodo,
        "C": methane,
        "tBu": tert_butyl,
        "phenyl": phenyl,
    }
    for atom_id, substituent in (substituents or {}).items():
        if substituent not in guests:
            raise ValueError(f"Unsupported BINOL substituent: {substituent}")
        result = _replace_with_fragment(
            result, atom_id, guests[substituent], name=name
        )
    result.name = name
    return result


binol = _binol("BINOL")
binol_3_3_dibromo = _binol("3,3'-dibromo-BINOL", {2: "Br", 13: "Br"})
binol_3_3_diiodo = _binol("3,3'-diiodo-BINOL", {2: "I", 13: "I"})
binol_6_6_dibromo = _binol("6,6'-dibromo-BINOL", {8: "Br", 19: "Br"})
binol_3_3_dimethyl = _binol("3,3'-dimethyl-BINOL", {2: "C", 13: "C"})


def _compose_binaphthyl_base() -> Structure:
    """Join two naphthalene fragments through their 1-position ports."""
    join_port = _port_at_atom(naphthalene, 1)
    return substitute(
        naphthalene,
        {join_port: (naphthalene, join_port)},
        name="1,1'-binaphthyl",
    )


def _binol_diformyl() -> Structure:
    """Build 2,2'-diformyl-1,1'-binaphthyl from naphthalene and carbonyl fragments."""
    result = _compose_binaphthyl_base()
    result = _replace_with_fragment(
        result, 2, carbonyl, name="2-formyl-1,1'-binaphthyl"
    )
    result = _replace_with_fragment(
        result, 13, carbonyl, name="2,2'-diformyl-1,1'-binaphthyl"
    )
    return result


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
        if amine == "2-aminopyridine":
            edges.extend([
                (imine_n, start + 1, BondType.SINGLE),
                (start, start + 1, BondType.SINGLE),
                (start + 1, start + 2, BondType.DOUBLE),
                (start + 2, start + 3, BondType.SINGLE),
                (start + 3, start + 4, BondType.DOUBLE),
                (start + 4, start + 5, BondType.SINGLE),
                (start + 5, start, BondType.DOUBLE),
            ])
        else:
            # 2-aminopyrimidine: attach at C2 and include N at ring position 3.
            edges.extend([
                (imine_n, start + 1, BondType.SINGLE),
                (start, start + 1, BondType.DOUBLE),
                (start + 1, start + 2, BondType.SINGLE),
                (start + 2, start + 3, BondType.DOUBLE),
                (start + 3, start + 4, BondType.SINGLE),
                (start + 4, start + 5, BondType.DOUBLE),
                (start + 5, start, BondType.SINGLE),
            ])
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
