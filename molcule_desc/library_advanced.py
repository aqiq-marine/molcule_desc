"""Higher-level reusable fragments composed from basic library structures."""

from .core import Port, PortRef, Structure, bond, ring, single, substitute
from .library import benzene, hydroxyl, methane, naphthol, phenyl, propane


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


def _phenol() -> Structure:
    """Build phenol from the reusable benzene and hydroxyl fragments.

    Ports 4 and 5 are kept on adjacent ring atoms so a four-carbon chain can
    be attached across them and closed into the second ring.
    """
    return substitute(
        benzene,
        {single(3): (hydroxyl, single(1))},
        expose={
            single(1): (0, single(1)),
            single(2): (0, single(2)),
            single(3): (0, single(4)),
            single(4): (0, single(5)),
            single(5): (0, single(6)),
        },
        name="phenol",
    )


phenol = _phenol()


def _butane_chain() -> Structure:
    """Build a four-carbon chain from propane and methane with both ends open."""
    return bond(
        propane,
        single(1),
        methane,
        single(1),
        expose={
            single(1): PortRef(0, single(6)),
            single(2): PortRef(1, single(2)),
        },
        name="butane chain",
    )


def _tetralol() -> Structure:
    """Build tetralol by attaching a butane chain to phenol and closing the ring."""
    chain = _butane_chain()
    attached = bond(
        phenol,
        single(4),
        chain,
        single(1),
        expose={
            single(1): PortRef(0, single(5)),
            single(2): PortRef(1, single(2)),
        },
        name="phenol with butane chain",
    )
    return ring(
        attached,
        single(1),
        single(2),
        expose={},
        name="tetralol",
    )

tetralol = _tetralol()


def _compose_h8_binol() -> Structure:
    """Join two tetralol fragments at their 1-positions (H8-BINOL)."""
    join_port = _port_at_atom(tetralol, 1)
    return substitute(
        tetralol, {join_port: (tetralol, join_port)}, name="H8-BINOL"
    )


h8_binol = _compose_h8_binol()
