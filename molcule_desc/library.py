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
    # Kekule representation: alternating double and single bonds around a six-membered ring.
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
        ports={single(1): 0, single(2): 5},
    )
    return ring(chain, single(1), single(2), expose={}, name="benzene")


benzene = _benzene()
