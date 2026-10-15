"""Higher-level reusable fragments composed from basic library structures."""

from .core import BondType, Port, Structure, single, substitute
from .library import _graph_structure, methane, naphthol, phenyl


def _tert_butyl() -> Structure:
    """Build tert-butyl from methane fragments, exposing one attachment port."""
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
    return substitute(phenyl, {single(1): (phenyl, single(1))}, name="biphenyl")


biphenyl = _biphenyl()


def _port_at_atom(structure: Structure, atom_id: int) -> Port:
    """Return an exposed port anchored to a specified atom."""
    for port, anchor in structure.ports.items():
        if anchor == atom_id:
            return port
    raise ValueError(f"{structure.name} has no exposed port on atom {atom_id}")


def _compose_binol_base() -> Structure:
    """Compose BINOL from two naphthol fragments joined at their 1-positions."""
    join_port = _port_at_atom(naphthol, 1)
    return substitute(naphthol, {join_port: (naphthol, join_port)}, name="BINOL")


binol = _compose_binol_base()


def _tetralol() -> Structure:
    """Build a 2-tetralol fragment with the phenolic ring remaining aromatic.

    Atom numbering follows the naphthol fragment: the hydroxyl-bearing
    aromatic ring is atoms 0-5; atoms 4-5 is the shared fused bond. The
    phenolic ring keeps its aromatic Kekule pattern. The opposite ring is
    hydrogenated by changing its non-fused bonds to single bonds.
    """
    return _graph_structure(
        "tetralol",
        ["C"] * 10 + ["O"],
        [
            (0, 1, BondType.DOUBLE), (1, 2, BondType.SINGLE),
            (2, 3, BondType.DOUBLE), (3, 4, BondType.SINGLE),
            (4, 5, BondType.SINGLE), (5, 0, BondType.DOUBLE),
            (5, 6, BondType.SINGLE), (6, 7, BondType.SINGLE),
            (7, 8, BondType.SINGLE), (8, 9, BondType.SINGLE),
            (9, 4, BondType.SINGLE), (2, 10, BondType.SINGLE),
        ],
        valences=[4] * 10 + [2],
    )


tetralol = _tetralol()


def _compose_h8_binol() -> Structure:
    """Join two tetralol fragments at their 1-positions (H8-BINOL)."""
    join_port = _port_at_atom(tetralol, 1)
    return substitute(
        tetralol, {join_port: (tetralol, join_port)}, name="H8-BINOL"
    )


h8_binol = _compose_h8_binol()
