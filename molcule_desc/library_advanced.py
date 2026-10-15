"""Higher-level reusable fragments composed from basic library structures."""

from .core import Atom, Bond, BondType, Port, PortRef, Structure, bond, single
from .library import methane


def _tert_butyl() -> Structure:
    """Build tert-butyl by attaching three methane fragments to one methane.

    One port on the central carbon is retained as the substituent attachment;
    the remaining ports on the three terminal carbons become their hydrogens
    when the completed molecule is requested.
    """
    central = methane
    next_port = 1
    for i in range(3):
        methyl = methane
        # The central carbon's first three ports are used to attach methyls.
        central_port = single(next_port)
        methyl_port = single(1)
        remaining: dict[Port, PortRef] = {}
        out_index = 1
        for port in central.ports:
            if port != central_port:
                remaining[single(out_index)] = PortRef(0, port)
                out_index += 1
        for port in methyl.ports:
            if port != methyl_port:
                remaining[single(out_index)] = PortRef(1, port)
                out_index += 1
        central = bond(
            central, central_port, methyl, methyl_port,
            expose=remaining, name=f"tert-butyl assembly {i + 1}",
        )
        next_port = 1
    # The first remaining port is the central carbon's fourth valence.
    # All other ports are the three hydrogens on each methyl carbon.
    # Keep only the central carbon port exposed; methyl ports are internal
    # terminal valences and will be hydrogen-completed on finalization.
    central_atom_id = next(
        atom_id for atom_id, atom in central.atoms.items()
        if atom.element == "C"
        and sum(1 for edge in central.bonds if atom_id in (edge.start, edge.end)) == 3
    )
    attachment_port = next(
        port for port, atom_id in central.ports.items() if atom_id == central_atom_id
    )
    atoms = dict(central.atoms)
    bonds = list(central.bonds)
    next_atom_id = max(atoms) + 1
    for port, anchor in central.ports.items():
        if anchor == central_atom_id:
            continue
        hydrogen_id = next_atom_id
        next_atom_id += 1
        atoms[hydrogen_id] = Atom(hydrogen_id, "H", 1)
        bonds.append(Bond(anchor, hydrogen_id, BondType.SINGLE))
    return Structure(
        "tert-butyl",
        atoms.values(),
        bonds,
        {single(1): central_atom_id},
    )




tert_butyl = _tert_butyl()
