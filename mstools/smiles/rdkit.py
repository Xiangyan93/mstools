#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Dict, Iterator, List, Optional, Union, Literal, Tuple
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions
)


def smiles2inchi(smiles: str) -> str:
    rdk_mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchi(rdk_mol)


def inchi2smiles(inchi: str) -> str:
    rdk_mol = Chem.MolFromInchi(inchi)
    return Chem.MolToSmiles(rdk_mol)


def get_charge(smiles: str) -> int:
    rdk_mol = Chem.MolFromSmiles(smiles)
    return Chem.GetFormalCharge(rdk_mol)


def get_heavy_atoms_number(smiles: str) -> int:
    rdk_mol = Chem.MolFromSmiles(smiles)
    return rdk_mol.GetNumAtoms()


def get_rdkit_smiles(smiles):
    rdk_mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(rdk_mol)


def get_stereo_isomer(smiles: str) -> List[str]:
    rdkit_mol = Chem.MolFromSmiles(smiles)
    isomers = tuple(EnumerateStereoisomers(rdkit_mol))
    return list(map(Chem.MolToInchi, isomers))


def add_atom(smiles: str,
             add_atom_number: int,
             add_bond_order: int,
             atom_number: str = 6) -> List[str]:
    """ Generate a set of molecules by adding an atom bonded to input molecule.

    Parameters
    ----------
    smiles: str
        SMILES string of input molecule.
    add_atom_number: int
        The atomic number of the atom to be added.
    add_bond_order: int
        The bond order of the bond to be added.
    atom_number: int
        The atomic number of the atoms of the input molecules that are possible to connect to the added atom.
        0 means use all atoms will be tried.

    Returns
    -------
    A list of SMILES strings.
    """
    new_smiles_list = []
    rdk_mol = Chem.MolFromSmiles(smiles)
    for atom in rdk_mol.GetAtoms():
        if atom.GetAtomicNum() != atom_number:
            continue
        if atom.GetImplicitValence() < add_bond_order:
            continue

        new_mol = Chem.RWMol(rdk_mol)
        new_id = new_mol.AddAtom(Chem.Atom(add_atom_number))
        if add_bond_order == 1:
            new_mol.AddBond(atom.GetIdx(), new_id, Chem.BondType.SINGLE)
        elif add_bond_order == 2:
            new_mol.AddBond(atom.GetIdx(), new_id, Chem.BondType.DOUBLE)
        elif add_bond_order == 3:
            new_mol.AddBond(atom.GetIdx(), new_id, Chem.BondType.TRIPLE)
        else:
            continue
        new_smiles = get_rdkit_smiles(Chem.MolToSmiles(new_mol))
        if new_smiles not in new_smiles_list:
            new_smiles_list.append(new_smiles)
    return new_smiles_list


def replace_atom(smiles: str,
                 old_atom_number: int,
                 new_atom_number: int) -> List[str]:
    """ Generate a set of molecules by replacing an atom by another atom in input molecule.

    Parameters
    ----------
    smiles: str
        SMILES string of input molecule.
    old_atom_number: int
        The atomic number of the atom to be replaced.
    new_atom_number
        The atomic number of the replaced atom.
    Returns
    -------
    A list of SMILES strings.
    """
    new_smiles_list = []
    rdk_mol = Chem.MolFromSmiles(smiles)
    for atom in rdk_mol.GetAtoms():
        if atom.GetAtomicNum() != old_atom_number:
            continue

        new_mol = Chem.RWMol(rdk_mol)
        idx = atom.GetIdx()
        new_mol.GetAtoms()[idx].SetAtomicNum(new_atom_number)
        try:
            new_smiles = get_rdkit_smiles(Chem.MolToSmiles(new_mol))
        except:
            continue
        else:
            if new_smiles not in new_smiles_list:
                new_smiles_list.append(new_smiles)
    return new_smiles_list
