#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Dict, Iterator, List, Optional, Union, Literal, Tuple
import os


class Packmol:
    """Python API for Packmol. It is used to build initial box for molecular simulation.

    Parameters
    ----------
    exe: str
        executable file of packmol.
    """
    def __init__(self, exe: str):
        self.exe = exe

    def build_box(self, pdb_files: List[str], n_mol_list: List[int], output_file: str,
                  box_size: Tuple[float, float, float],
                  slab: float = None,
                  tolerance: float = 1.8, seed: int = 0,
                  inp_file='build.inp', silent: bool = False):
        """

        Parameters
        ----------
        pdb_files: List[str]
            A list of pdb files, each of them contains only single molecule.
        n_mol_list: List[int]
            The number of molecules to be created.
        output_file: str
            The name of output file.
        box_size: Tuple[float, float, float]
            The length of simulation box, 3 floats in Angstrom.
        slab: float or None
            If float, a vacuum space in Z direction will be created.
            And the depth of the vacuum space is assigned by slab.
        tolerance: float
            Tolerance of packmol.
        seed: int
            Random seed.
        inp_file: str
            The name of packmol input file.
        silent: bool
            If False, more details will presented in screen.
        Returns
        -------

        """
        if slab is None:
            inp = self.__build_bulk(pdb_files=pdb_files,
                                    n_mol_list=n_mol_list,
                                    output=output_file,
                                    box_size=tuple([l - 1.0 for l in box_size]),
                                    tolerance=tolerance,
                                    seed=seed)
        else:
            inp = self.__build_slab(pdb_files=pdb_files,
                                    n_mol_list=n_mol_list,
                                    output=output_file,
                                    box_size=[l - 1.0 for l in box_size],
                                    slab=slab,
                                    tolerance=tolerance,
                                    seed=seed)

        with open(inp_file, 'w') as f:
            f.write(inp)

        # TODO subprocess PIPE not work for Packmol new version, do not know why.
        if silent:
            if os.name == 'nt':
                os.system(self.exe + ' < %s > nul' % inp_file)
            else:
                os.system(self.exe + ' < %s > /dev/null' % inp_file)
        else:
            os.system(self.exe + ' < %s' % inp_file)

            # (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            # sp = subprocess.Popen([self.PACKMOL_BIN], stdin=PIPE, stdout=stdout, stderr=stderr)
            # sp.communicate(input=inp.encode())

    def __build(self, pdb_files: List[str], n_mol_list: List[int],
                output: str, tolerance: float = 1.8, seed: int = 0):
        if len(pdb_files) == 0:
            raise RuntimeError(f'Empty list of pdb_files as input. {pdb_files}')
        if len(pdb_files) != len(n_mol_list):
            raise RuntimeError(f'pdb_files and n_mol_list must have the same length. '
                               f'{len(pdb_files)}, {len(n_mol_list)}.')

        filetypes = {filename.split('.')[-1].lower() for filename in pdb_files}
        if len(filetypes) > 1:
            raise RuntimeError(f'All file types must be the same. {pdb_files}')
        filetype = filetypes.pop()

        inp = (
            'filetype {filetype}\n'
            'tolerance {tolerance}\n'
            'output {output}\n'
            'seed {seed}\n'
            'add_box_sides 1.0\n'.format(filetype=filetype, tolerance=tolerance, output=output, seed=seed)
        )
        return inp

    def __build_bulk(self, pdb_files: List[str], n_mol_list: List[int], output: str,
                     box_size: Tuple[float, float, float],
                     tolerance: float = 1.8, seed: int = 0) -> str:
        inp = self.__build(pdb_files=pdb_files, n_mol_list=n_mol_list, output=output, tolerance=tolerance, seed=seed)
        for i, filename in enumerate(pdb_files):
            number = n_mol_list[i]
            box = '0 0 0 %f %f %f' % tuple(box_size)
            inp += (
                'structure {filename}\n'
                'number {number}\n'
                'inside box {box}\n'
                'end structure\n'.format(filename=filename, number=number, box=box)
            )
        return inp

    def __build_slab(self, pdb_files: List[str], n_mol_list: List[int], output: str,
                     box_size: Tuple[float, float, float], slab: float,
                     tolerance: float = 1.8, seed: int = 0) -> str:
        inp = self.__build(pdb_files=pdb_files, n_mol_list=n_mol_list, output=output, tolerance=tolerance, seed=seed)
        box_liq = '0 0 0 %f %f %f' % (box_size[0], box_size[1], slab)
        box_gas = '0 0 %f %f %f %f' % (slab, box_size[0], box_size[1], box_size[2])
        for i, filename in enumerate(pdb_files):
            # put 1/50 molecules in gas phase. Do not put too many in case of nucleation in gas phase
            n_gas = n_mol_list[i] // 50
            n_liq = n_mol_list[i] - n_gas
            inp += (
                'structure {filename}\n'
                'number {n_liq}\n'
                'inside box {box_liq}\n'
                'end structure\n'
                'structure {filename}\n'
                'number {n_gas}\n'
                'inside box {box_gas}\n'
                'end structure\n'.format(filename=filename, n_liq=n_liq, n_gas=n_gas,
                                         box_liq=box_liq, box_gas=box_gas)
            )
        return inp
