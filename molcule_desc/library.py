"""Reusable elementary structures and example compounds."""

from .core import (
    Atom,
    Bond,
    BondType,
    Port,
    PortRef,
    Structure,
    bond,
    double,
    single,
    triple,
)

# Terminal/dummy structures. XS is a terminal H; XD is a CH2-side fragment;
# XT is a CH-side fragment. XD/XT expose their hydrogen positions as ports.
XS = Structure.atom("H", [single(1)], name="XS")
XD = Structure.atom("C", [double(1), single(1), single(2)], name="XD")
XT = Structure.atom("C", [triple(1), single(1)], name="XT")


def _ethylene() -> Structure:
    return bond(
        XD, double(1), XD, double(1),
        expose={
            single(1): PortRef(0, single(1)),
            single(2): PortRef(0, single(2)),
            single(3): PortRef(1, single(1)),
            single(4): PortRef(1, single(2)),
        },
        name="ethylene",
    )


def _methanimine() -> Structure:
    carbon = Structure.atom("C", [double(1), single(1), single(2)], name="imine C")
    nitrogen = Structure.atom("N", [double(1), single(1)], name="imine N")
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
    carbon = Structure.atom("C", [double(1), single(1), single(2)], name="carbonyl C")
    oxygen = Structure.atom("O", [double(1)], name="carbonyl O")
    return bond(
        carbon, double(1), oxygen, double(1),
        expose={
            single(1): PortRef(0, single(1)),
            single(2): PortRef(0, single(2)),
        },
        name="carbonyl",
    )


def _acetone() -> Structure:
    carbonyl = _carbonyl()
    methyl = Structure.atom(
        "C", [single(1), single(2), single(3), single(4)], name="methyl C"
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
acetone = _acetone()
