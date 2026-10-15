""""Higher-level reusable fragments composed from basic library structures."""

from .core import Port, Structure, single, substitute
from .library import methane, naphthol, phenyl


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
    """Compose BINOL from two naphthol fragments joined at their 1-positions.

    The two consumed biaryl ports are the only ports removed. Every other
    naphthol port, including ring-substitution positions and both O-H ports,
    remains exposed so that this scaffold can be reused compositionally.
    """
    join_port = _port_at_atom(naphthol, 1)
    return substitute(naphthol, {join_port: (naphthol, join_port)}, name="BINOL")


binol = _compose_binol_base()
