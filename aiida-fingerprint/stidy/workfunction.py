# -*- coding: utf-8 -*-

import re
from os import path, remove
from string import ascii_uppercase
from distutils.spawn import find_executable
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE, STDOUT
from pymatgen import Structure, Lattice
from pymatgen.io.cif import CifWriter
from aiida.orm import DataFactory
from aiida.work.workfunctions import workfunctions

StructureData = DataFactory('structure')
CifData = DataFactory('cif')


@workfunctions
def stidy_workfunction(structure, ang=5, d1=0.55, d2=0.55, d3=0.55):
    def stidy(structure, ang, d1, d2, d3):
        PLATON = find_executable('platon')
        if not PLATON:
            PLATON = '../bin/platon'

        with NamedTemporaryFile(suffix='.cif') as temp_file:
            # write temporary cif file
            CifWriter(structure).write_file(temp_file.name)
            temp_file.flush()
            # run ADDSYM_SHX to make PLATON recognize symmetries
            addsym_shx_process = Popen(['platon', '-o', temp_file.name],
                                       stdout=PIPE,
                                       stderr=STDOUT,
                                       stdin=PIPE)
            addsym_shx_process.communicate(
                input=b'ADDSYM_SHX {} {} {} {}'.format(ang, d1, d2, d3))
            # call STIDY on the ADDSYM_SHX output
            temp_file_dirname, temp_file_basena dme = path.split(
                temp_file.name)
            temp_file_basename_extless, _ = path.splitext(
                temp_file_basename)
            temp_file_basename_spf = temp_file_basename_extless + '_pl.spf'
            temp_file_spf = path.join(
                temp_file_dirname, temp_file_basename_spf)

            stidy_process = Popen(['platon', '-o', temp_file_spf],
                                  stdout=PIPE,
                                  stderr=STDOUT,
                                  stdin=PIPE)
            stidy_data = stidy_process.communicate(input=b'STIDY')
        stidy_output = stidy_data[0].decode('utf-8')

        # clean up files
        if path.isfile('check.def'):
            remove('check.def')

        return stidy_output

    def get_space_group(stidy_output):
        regexp = re.compile(r'\s*Number in IT :\s*(\d+)')
        match = regexp.search(stidy_output)
        space_group = int(match.group(1))
        return space_group

    def get_wyckoffs(stidy_output):
        wyckoffs = []
        for site in self.ctx.sites:
            wyckoffs.append(site[1].split('(')[-1].split(')')[0])
        return wyckoffs

    def get_sites(stidy_output:
        regexp = re.compile(
            r'\s+([a-zA-Z]{1,2})(\d)+\s+([\w\(\)]{4,5})\s+([\d\/\.]+)\s+([\d\/\.]+)\s+([\d\/\.]+)\s+(\w+)\s+(\d+)')
        match_blocks = []
        for block in stidy_output.split('Wyckoff'):
            matches = []
            matches = regexp.findall(block)
            for m, match in enumerate(matches):
                match = list(match)
                for i in range(3, 6):
                    if '/' in match[i]:
                        j = [float(k) for k in match[i].split('/')]
                        num = j[0] / j[1]
                        match[i] = num
                    else:
                        match[i] = float(match[i])
                match[0] = ''.join(match[0:2])
                del match[1]
                match[6] = int(match[6])
                matches[m] = tuple(match)
            if matches:
                match_blocks.append(matches)
        sites = match_blocks[0]
        return sites

    def get_fingerprint(space_group, wyckoffs, sites):
        ascii_iter = iter(ascii_uppercase)
        symbols = [site[5] for site in sites]
        symbol_map = {}
        for symbol in symbols:
            if symbol not in symbol_map.keys():
                symbol_map[symbol] = next(ascii_iter)
        sites = [(symbol, wyckoff)
                 for symbol, wyckoff in zip(symbols, wyckoffs)]
        formatted_wyckoffs = []
        for wyckoff in sorted(set(wyckoffs)):
            formatted_wyckoff = []
            for symbol in sorted(set(symbols)):
                count = sites.count((symbol, wyckoff))
                if count:
                    formatted_wyckoff.append(
                        '{:02d}{}'.format(count, symbol_map[symbol]))
            formatted_wyckoffs.append(''.join(formatted_wyckoff + [wyckoff]))
        fingerprint = '_'.join([str(space_group)] + formatted_wyckoffs)
        return fingerprint

    ############################################################################

    if isinstance(structure, StructureData):
        structure = structure.get_pymatgen()
    elif isinstance(structure, Atoms):
        structure = AseAtomsAdaptor.get_structure(structure)
    # elif isinstance(structure, CifData):
        # structure = CifParser(CifData.get_attr(
        #     'source')).get_structures(primitive=False)[0]

    stidy_output = stidy(structure, ang, d1, d2, d3)
    space_group = get_space_group(stidy_output)
    wyckoffs = get_wyckoffs(stidy_output)
    sites = get_sites(stidy_output)
    fingerprint = get_fingerprint(space_group, wyckoffs, sites)
    return space_group, wyckoffs, sites, fingerprint
