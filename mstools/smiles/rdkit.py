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
