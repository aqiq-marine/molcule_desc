"""RDKit conversion utilities for molcule_desc structures.

RDKit is an optional dependency. Import it lazily so the core DSL remains
usable without RDKit installed.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit.Chem.rdchem import Mol

from .core import BondType, Structure


def to_rdkit_mol(
    structure: Structure,
    *,
    finalize: bool = True,
) -> "Mol":
    """Convert a Structure to an RDKit Mol.

    By default, exposed ports are completed using Structure.finalize() before
    conversion. Set finalize=False to require a structure with no exposed ports.
    The returned RDKit molecule is sanitized, but stereochemistry and aromatic
    bond types are not inferred by this converter.
    """
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError(
            "RDKit is required for molecular conversion. "
            "Install it with pip install 'molcule-desc[chemistry]'."
        ) from exc

    if finalize:
        structure = structure.finalize()
    elif structure.ports:
        raise ValueError(
            "Structure has exposed ports; call finalize() or pass finalize=True."
        )

    rw_mol = Chem.RWMol()
    atom_indices: dict[int, int] = {}

    for atom_id, atom in structure.atoms.items():
        rd_atom = Chem.Atom(atom.element)
        rd_atom.SetFormalCharge(atom.formal_charge)
        atom_indices[atom_id] = rw_mol.AddAtom(rd_atom)

    bond_types = {
        BondType.SINGLE: Chem.BondType.SINGLE,
        BondType.DOUBLE: Chem.BondType.DOUBLE,
        BondType.TRIPLE: Chem.BondType.TRIPLE,
    }

    for bond in structure.bonds:
        rw_mol.AddBond(
            atom_indices[bond.start],
            atom_indices[bond.end],
            bond_types[bond.order],
        )

    mol = rw_mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def to_smiles(
    structure: Structure,
    *,
    canonical: bool = True,
    isomeric: bool = True,
    finalize: bool = True,
) -> str:
    """Convert a Structure to a SMILES string using RDKit.

    Explicit hydrogen atoms used internally by molcule_desc are removed where
    RDKit considers them removable, producing conventional compact SMILES.
    """
    mol = to_rdkit_mol(structure, finalize=finalize)
    from rdkit import Chem

    mol = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(
        mol,
        canonical=canonical,
        isomericSmiles=isomeric,
    )
