"""Examples of lactone formation and iodoetherification using bond() then ring()."""
from molcule_desc.core import PortRef, atom, bond, ported_atom, ring, single
from molcule_desc.library import carbonyl, iodo


def _carbon(name):
    return ported_atom(atom("C", 4), [single(i) for i in range(1, 5)], name=name)


def _oxygen(name):
    return ported_atom(atom("O", 2), [single(1), single(2)], name=name)


def _join(left, lp, right, rp, expose, name):
    return bond(left, lp, right, rp, expose=expose, name=name)


def _lactone(ring_size):
    """Build the open ester chain with bond(), then close it with ring()."""
    if not 4 <= ring_size <= 8:
        raise ValueError("Lactone ring size must be between 4 and 8")

    chain = _join(
        carbonyl, single(1), _carbon("methylene"), single(1),
        {single(1): PortRef(0, single(2)), single(2): PortRef(1, single(2))},
        "open lactone chain",
    )
    active = single(2)
    for index in range(ring_size - 3):
        chain = _join(
            chain, active, _carbon(f"methylene {index + 2}"), single(1),
            {single(1): PortRef(0, single(1)), single(2): PortRef(1, single(2))},
            "open lactone chain",
        )
        active = single(2)

    chain = _join(
        chain, active, _oxygen("lactone oxygen"), single(1),
        {single(1): PortRef(0, single(1)), single(2): PortRef(1, single(2))},
        "open lactone chain",
    )
    return ring(chain, single(1), single(2), expose={}, name=f"{ring_size}-membered lactone")


beta_propiolactone = _lactone(4)
gamma_butyrolactone = _lactone(5)
delta_valerolactone = _lactone(6)
epsilon_caprolactone = _lactone(7)
seven_membered_lactone = _lactone(8)


def _iodoetherification_product():
    """Build 2-(iodomethyl)tetrahydrofuran as a chain and close its ring."""
    chain = _join(
        _oxygen("ring oxygen"), single(1), _carbon("substituted ring carbon"), single(1),
        {
            single(1): PortRef(0, single(2)),
            single(2): PortRef(1, single(2)),
            single(3): PortRef(1, single(3)),
        },
        "open iodoether chain",
    )
    branch = _join(
        _carbon("iodomethyl carbon"), single(1), iodo, single(1),
        {single(1): PortRef(0, single(2))},
        "iodomethyl branch",
    )
    chain = _join(
        chain, single(3), branch, single(1),
        {single(1): PortRef(0, single(1)), single(2): PortRef(0, single(2))},
        "open iodoether chain",
    )
    for index in range(3):
        chain = _join(
            chain, single(2), _carbon(f"ring methylene {index + 1}"), single(1),
            {single(1): PortRef(0, single(1)), single(2): PortRef(1, single(2))},
            "open iodoether chain",
        )
    return ring(
        chain, single(1), single(2), expose={},
        name="2-(iodomethyl)tetrahydrofuran",
    )


iodoetherification_5_exo = _iodoetherification_product()
