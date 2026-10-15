"""Ring-forming reaction examples built with the core graph API.

The lactones are cyclic ester products with ring sizes counted by atoms in
the ring (including the ester oxygen and carbonyl carbon). The iodoether is a
representative 5-exo cyclization product, 2-(iodomethyl)tetrahydrofuran.

These are structure-construction examples; they do not model reaction
conditions, regiochemical selectivity, or a reaction mechanism.
"""
from __future__ import annotations

from molcule_desc.core import Atom, Bond, BondType, Structure, ring, single


def _lactone(ring_size: int) -> Structure:
    """Build an n-membered saturated lactone by closing a linear precursor."""
    if not 4 <= ring_size <= 8:
        raise ValueError("Lactone ring size must be between 4 and 8")

    # Ring path: carbonyl C - (ring_size - 2) saturated carbons - O.
    # Atom 1 is the exocyclic carbonyl oxygen.
    atoms = [
        Atom(0, "C", 4),
        Atom(1, "O", 2),
        *[Atom(i, "C", 4) for i in range(2, ring_size)],
        Atom(ring_size, "O", 2),
    ]
    bonds = [
        Bond(0, 1, BondType.DOUBLE),
        Bond(0, 2, BondType.SINGLE),
    ]
    bonds.extend(
        Bond(i, i + 1, BondType.SINGLE)
        for i in range(2, ring_size)
    )
    # The carbonyl carbon (0) connects to the first ring carbon (2); the
    # chain continues through the methylene units to ring oxygen.
    ports = {single(1): 0, single(2): ring_size}
    precursor = Structure(
        f"{ring_size}-membered lactone precursor", atoms, bonds, ports
    )
    return ring(
        precursor,
        single(1),
        single(2),
        expose={},
        name=f"{ring_size}-membered lactone",
    )


# Four-membered beta-propiolactone through eight-membered lactones.
beta_propiolactone = _lactone(4)
gamma_butyrolactone = _lactone(5)
delta_valerolactone = _lactone(6)
epsilon_caprolactone = _lactone(7)
seven_membered_lactone = _lactone(8)


def _iodoetherification_product() -> Structure:
    """Build 2-(iodomethyl)tetrahydrofuran as a 5-exo iodoetherification product.

    The ring is O-C-C-C-C. Iodine is attached to the exocyclic methylene
    carbon that corresponds to the terminal alkene carbon in a 4-pentenol
    cyclization.
    """
    atoms = [
        Atom(0, "O", 2),  # ring oxygen
        Atom(1, "C", 4),  # substituted ring carbon
        Atom(2, "C", 4),
        Atom(3, "C", 4),
        Atom(4, "C", 4),  # terminal ring carbon
        Atom(5, "C", 4),  # exocyclic CH2I
        Atom(6, "I", 1),
    ]
    bonds = [
        Bond(0, 1, BondType.SINGLE),
        Bond(1, 2, BondType.SINGLE),
        Bond(2, 3, BondType.SINGLE),
        Bond(3, 4, BondType.SINGLE),
        Bond(1, 5, BondType.SINGLE),
        Bond(5, 6, BondType.SINGLE),
    ]
    precursor = Structure(
        "iodoetherification precursor",
        atoms,
        bonds,
        ports={single(1): 0, single(2): 4},
    )
    return ring(
        precursor,
        single(1),
        single(2),
        expose={},
        name="2-(iodomethyl)tetrahydrofuran",
    )


iodoetherification_5_exo = _iodoetherification_product()
