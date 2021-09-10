#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Dict, Iterator, List, Optional, Union, Literal, Tuple
import os
import shutil
from .packmol import Packmol
from ..molecule import Molecule


def build_cubic_box(packmol_exe: str, work_dir: str, file_type: Literal['pdb', 'xyz'],
                    smiles_list: List[str], n_mol_list: List[int], density: float,
                    save_tmp_file: bool = False):
    """This function is used to build a cubic box of molecules.

    Parameters
    ----------
    packmol_exe: str
        Executable file of packmol.
    work_dir: str
        The working directory.
    file_type: Literal['pdb', 'xyz']
        output file type.
    smiles_list: List[str]
        A list of SMILES ings.
    n_mol_list: List[int]
        How many molecules to be created. This is corresponding to smiles_list.
    density: float
        the density
    Returns
    -------

    """
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    cwd = os.getcwd()
    os.chdir(work_dir)
    pdb_files = []
    n_components = len(smiles_list)
    molwt_list = []  # molecule weight of each molecule
    # generate 3D-structure files of single molecules.
    for i, smiles in enumerate(smiles_list):
        pdb = 'mol-%i.%s' % (i, file_type)
        # mol2 = 'mol-%i.mol2' % i
        mol3d = Molecule(smiles)
        mol3d.write(pdb, filetype=file_type)
        # mol3d.write(mol2, filetype='mol2')
        pdb_files.append(pdb)
        molwt_list.append(mol3d.molwt)
    # calculate the box size based on density and molecular weights.
    mass = sum([molwt_list[i] * n_mol_list[i] for i in range(n_components)])
    vol = 10 / 6.022 * mass / density
    length = vol ** (1 / 3)
    box_size = (length, length, length)

    Packmol(exe=packmol_exe).build_box(
        pdb_files=pdb_files,
        n_mol_list=n_mol_list,
        output_file='packmol.%s' % file_type,
        box_size=box_size,
        silent=True)

    os.chdir(cwd)
    shutil.copy(os.path.join(work_dir, 'packmol.%s' % file_type), os.path.join(cwd))
    if not save_tmp_file:
        shutil.rmtree(work_dir)
