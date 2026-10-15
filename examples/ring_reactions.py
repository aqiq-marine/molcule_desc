"""Examples of lactone formation and iodoetherification using library fragments."""

from molcule_desc.core import PortRef, bond, ring, single
from molcule_desc.library import carbonyl, hydroxyl, iodo, methylene


def _join(left, lp, right, rp, expose, name):
    """Join fragments and let bond() cap every port not explicitly exposed."""
    return bond(left, lp, right, rp, expose=expose, name=name)


def _lactone(ring_size):
    """Build the open ester chain with bond(), then close it with ring()."""
    if not 4 <= ring_size <= 8:
        raise ValueError("Lactone ring size must be between 4 and 8")

    # Keep one carbonyl-carbon port for ring closure and one chain-growth port.
    chain = _join(
        carbonyl, single(1), methylene, single(1),
        {
            single(1): PortRef(0, single(2)),
            single(2): PortRef(1, single(2)),
        },
        "open lactone chain",
    )
    active = single(2)

    # The first methylene is already present; add enough to reach the desired
    # ring size. Unexposed methylene ports are capped with hydrogen by bond().
    for _ in range(ring_size - 4):
        chain = _join(
            chain, active, methylene, single(1),
            {
                single(1): PortRef(0, single(1)),
                single(2): PortRef(1, single(2)),
            },
            "open lactone chain",
        )
        active = single(2)

    # Connect the terminal methylene to oxygen, preserving the oxygen's other
    # port and the carbonyl closure port. ring() then joins those two ports.
    chain = _join(
        chain, active, hydroxyl, single(1),
        {
            single(1): PortRef(0, single(1)),
            single(2): PortRef(1, single(2)),
        },
        "open lactone chain",
    )
    return ring(
        chain, single(1), single(2), expose={},
        name=f"{ring_size}-membered lactone",
    )


beta_propiolactone = _lactone(4)
gamma_butyrolactone = _lactone(5)
delta_valerolactone = _lactone(6)
epsilon_caprolactone = _lactone(7)
seven_membered_lactone = _lactone(8)


def _iodoetherification_product():
    """Build 2-(iodomethyl)tetrahydrofuran as a chain and close its ring."""
    # Ring oxygen -> substituted ring carbon. Preserve the oxygen closure port,
    # the carbon's ring-growth port, and its substituent port.
    chain = _join(
        hydroxyl, single(1), methylene, single(1),
        {
            single(1): PortRef(0, single(2)),
            single(2): PortRef(1, single(2)),
            single(3): PortRef(1, single(3)),
        },
        "open iodoether chain",
    )

    # Attach the iodomethyl substituent; only the methylene's branch port is
    # exported, while all other unused ports are capped at this step.
    branch = _join(
        methylene, single(1), iodo, single(1),
        {single(1): PortRef(0, single(2))},
        "iodomethyl branch",
    )
    chain = _join(
        chain, single(3), branch, single(1),
        {
            single(1): PortRef(0, single(1)),
            single(2): PortRef(0, single(2)),
        },
        "open iodoether chain",
    )

    # Add three more ring methylenes and close the last one to the oxygen.
    for _ in range(3):
        chain = _join(
            chain, single(2), methylene, single(1),
            {
                single(1): PortRef(0, single(1)),
                single(2): PortRef(1, single(2)),
            },
            "open iodoether chain",
        )
    return ring(
        chain, single(1), single(2), expose={},
        name="2-(iodomethyl)tetrahydrofuran",
    )


iodoetherification_5_exo = _iodoetherification_product()
