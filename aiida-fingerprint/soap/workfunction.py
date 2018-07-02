# -*- coding: utf-8 -*-

import numpy as np
from sys import path
from aiida.orm import DataFactory

path.insert(0, 'GLOSIM2')
from libmatch.soap import get_soap
from libmatch.utils import ase2qp, get_spkit, get_spkitMax

StructureData = DataFactory('structure')
CifData = DataFactory('cif')


@workfunction
def soap_workfunction(aiida_structure, spkit_max, do_anonimization=True,
                      do_scaling=True, scale_per='site', **soapargs):
    def anonymize(structure):
        n_atoms = structure.n
        structure.set_atomic_numbers([1] * n_atoms)
        structure.set_chemical_symbols(['H'] * n_atoms)
        return structure

    def scale(structure, per='site'):
        if per == 'site':
            n_atoms = structure.n
        elif per == 'cell':
            n_atoms = 1  # scaling volume to 1/cell is the same as having 1 atom
        else:
            # raise error here
            pass
        new_cell = structure.get_cell() / np.cbrt(structure.cell_volume() / n_atoms)
        structure.set_cell(new_cell)
        new_pos = structure.get_positions() / \
            np.linalg.norm(structure.get_cell(), axis=1) * \
            np.linalg.norm(new_cell, axis=1)
        structure.set_positions(new_pos)
        return structure

    ############################################################################

    ase_atoms = aiida_structure.get_ase()
    quippy_atoms = ase2qp(ase_atoms)
    if do_anonimization:
        quippy_atoms = anonymize(quippy_atoms)
    if do_scaling:
        quippy_atoms = scale(quippy_atoms, per=scale_per)
    soap = get_soap(quippy_atoms, spkit=get_spkit(
        quippy_atoms), spkitMax=spkit_max, **soapargs)
    return soap
