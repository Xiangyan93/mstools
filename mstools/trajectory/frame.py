#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Dict, Iterator, List, Optional, Union, Literal, Tuple
from mendeleev import element


class PDB:
    def __init__(self, file: str):
        self.atomic_number = []
        self.positions = []
        for line in open(file, 'r').readlines():
            if line.startswith('CRYST1'):
                self.box_size = list(map(float, line.split()[1:4]))
                for degree in line.split()[4:6]:
                    assert degree == '90.00'
            elif line.startswith('HETATM'):
                info = line.split()
                self.atomic_number.append(element(info[-1]).atomic_number)
                self.positions.append(list(map(float, info[-6:-3])))
