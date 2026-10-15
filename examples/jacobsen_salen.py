"""Build the Jacobsen salen ligand shown in the project example image.

This example describes connectivity only. The current Structure model does not
encode tetrahedral stereochemistry, so the (R,R) configuration is not encoded.
Requires the optional RDKit integration for SMILES conversion.
"""

from molcule_desc import Atom, Bond, BondType, Structure, to_smiles


def build_jacobsen_salen() -> Structure:
    """Return N,N'-bis(3,5-di-tert-butylsalicylidene)-1,2-cyclohexanediamine.

    The two salicylidene rings are attached to imine carbons; each ring has
    an ortho hydroxy group and tert-butyl groups at the 3- and 5-positions.
    """
    atoms: list[Atom] = []
    bonds: list[Bond] = []

    def add_atom(element: str, valence: int = 4) -> int:
        atom_id = len(atoms)
        atoms.append(Atom(atom_id, element, valence))
        return atom_id

    def add_bond(a: int, b: int, order: BondType = BondType.SINGLE) -> None:
        bonds.append(Bond(a, b, order))

    def add_tert_butyl(anchor: int) -> None:
        central = add_atom("C")
        add_bond(anchor, central)
        for _ in range(3):
            methyl = add_atom("C")
            add_bond(central, methyl)

    imine_nitrogens: list[int] = []

    # Two 3,5-di-tert-butylsalicylidene halves.
    for _ in range(2):
        ring_c = [add_atom("C") for _ in range(6)]
        for i in range(6):
            # Kekulé benzene; RDKit perceives aromaticity on conversion.
            order = BondType.DOUBLE if i % 2 == 0 else BondType.SINGLE
            add_bond(ring_c[i], ring_c[(i + 1) % 6], order)

        # Ring position 0: imine substituent; position 1: OH;
        # positions 2 and 4: tert-butyl; positions 3 and 5: aromatic H.
        imine_c = add_atom("C")
        imine_n = add_atom("N", valence=3)
        add_bond(ring_c[0], imine_c)
        add_bond(imine_c, imine_n, BondType.DOUBLE)
        imine_nitrogens.append(imine_n)

        oxygen = add_atom("O", valence=2)
        add_bond(ring_c[1], oxygen)

        add_tert_butyl(ring_c[2])
        add_tert_butyl(ring_c[4])

    # 1,2-disubstituted cyclohexane: the two imine nitrogens attach to
    # adjacent stereogenic ring carbons. Chirality is not represented yet.
    cyclohexyl = [add_atom("C") for _ in range(6)]
    for i in range(6):
        add_bond(cyclohexyl[i], cyclohexyl[(i + 1) % 6])
    add_bond(imine_nitrogens[0], cyclohexyl[0])
    add_bond(imine_nitrogens[1], cyclohexyl[1])

    # Complete the two phenolic OH groups and all remaining valences with H.
    # Structure.finalize() also handles aromatic C-H and tert-butyl methyl H.
    # The O atoms above have one remaining single-bond valence each.
    ports = {}
    next_port = 1
    used = [0] * len(atoms)
    for edge in bonds:
        used[edge.start] += edge.order.value
        used[edge.end] += edge.order.value
    for atom_id, atom in enumerate(atoms):
        for _ in range(atom.valence - used[atom_id]):
            from molcule_desc import single
            ports[single(next_port)] = atom_id
            next_port += 1

    return Structure("Jacobsen salen ligand (connectivity only)", atoms, bonds, ports)


if __name__ == "__main__":
    ligand = build_jacobsen_salen()
    completed = ligand.finalize()
    print("Formula:", completed.formula_string())
    print("SMILES:", to_smiles(ligand))
