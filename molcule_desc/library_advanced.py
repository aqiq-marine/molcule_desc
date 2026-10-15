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
