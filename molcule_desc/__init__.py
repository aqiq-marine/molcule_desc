"""Compositional molecular structure prototype."""

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
    ring,
    substitute,
    double,
    single,
    triple,
)
from .smiles import to_rdkit_mol, to_smiles
from .library_advanced import biphenyl, tert_butyl
from .library import (
    hydrogen_single, carbon_double_two_single, carbon_triple_single, oxygen_double,
    acetone, benzene, carbonyl, ethylene, formaldehyde, methanimine,
    methane, ethane, propane, acetylene, cyclopropane, cyclobutane,
    cyclopentane, cyclohexane, cycloheptane, cyclooctane, cyclononane, cyclodecane, phenyl,
    hydroxyl, ether, carboxyl, ester, amino, amine, amide, nitrile, imine,
    pyridine, pyrrole, thiophene, imidazole, naphthalene,
    fluoro, chloro, bromo, iodo, thiol, thioether, sulfone, phosphine,
    phosphine_oxide, phosphate_ester, methylene, ethylene_linker, ethanol,
    diphenyl_ether, naphthol, indole, benzofuran, benzothiophene,
    quinoline, isoquinoline, anthracene,
)

__all__ = [
    "Atom", "Bond", "BondType", "Port", "PortRef", "Structure",
    "atom", "ported_atom", "bond", "ring", "substitute", "single", "double", "triple",
    "to_rdkit_mol", "to_smiles",
    "hydrogen_single", "carbon_double_two_single", "carbon_triple_single",
    "oxygen_double", "ethylene", "methanimine",
    "carbonyl", "formaldehyde", "acetone", "benzene", "tert_butyl",
    "methane", "ethane", "propane", "acetylene", "cyclopropane", "cyclobutane",
    "cyclopentane", "cyclohexane", "cycloheptane", "cyclooctane", "cyclononane", "cyclodecane",
    "phenyl", "hydroxyl", "ether", "carboxyl", "ester",
    "amino", "amine", "amide", "nitrile", "imine", "pyridine", "pyrrole",
    "thiophene", "imidazole", "naphthalene", "fluoro", "chloro", "bromo",
    "iodo", "thiol", "thioether", "sulfone", "phosphine",
    "phosphine_oxide", "phosphate_ester", "methylene", "ethylene_linker",
    "ethanol", "biphenyl", "diphenyl_ether", "naphthol", "indole",
    "benzofuran", "benzothiophene", "quinoline", "isoquinoline", "anthracene",
]
