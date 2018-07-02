# -*- coding: utf-8 -*-

import re
from os import path, remove
from string import ascii_uppercase
from distutils.spawn import find_executable
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE, STDOUT
from pymatgen import Structure, Lattice
from pymatgen.io.cif import CifWriter
from aiida.orm.data.base import Str, Float
from aiida.orm import DataFactory
from aiida.work.workchain import WorkChain
StructureData = DataFactory('structure')


class FingerprintWorkChainStidy(WorkChain):

    @classmethod
    def define(cls, spec):
        super(FingerprintWorkChainStidy, cls).define(spec)
        spec.input('structure', valid_type=StructureData)
        spec.outline(
            cls.get_pymatgen_structure,
            cls.stidy,
            cls.get_space_group,
            cls.get_sites,
            cls.get_wyckoffs,
            cls.get_cell,
            cls.get_stidy_structure,
            cls.get_fingerprint
        )
        spec.output('fingerprint', valid_type=Str)

    def get_pymatgen_structure(self):
        self.ctx.structure = StructureData.get_pymatgen(self.inputs.structure)

    def stidy(self):
        '''
        Run STRUCTURE TIDY as implemented in the PLATON software package.
        PLATON must either be in the PATH or in ../bin.

        References:
            A. L. Spek (2009). Acta Cryst., D65, 148-155.
            E. Parthé and L. M. Gelato (1984). Acta Cryst., A40, 169-183.

            L. M. Gelato and E. Parthé (1987). J. Appl. Cryst. 20, 139-143.
            S-Z. Hu and E. Parthé (2004). Chinese J. Struct. Chem. 23, 1150-1160.
        '''
        ang = 5
        d1 = 0.55
        d2 = 0.55
        d3 = 0.55

        PLATON = find_executable('platon')
        if not PLATON:
            PLATON = '../bin/platon'

        with NamedTemporaryFile(suffix='.cif') as temp_file:
            # write temporary cif file
            CifWriter(self.ctx.structure).write_file(temp_file.name)
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

        self.ctx.stidy_output = stidy_output

    def get_space_group(self):
        regexp = re.compile(r'\s*Number in IT :\s*(\d+)')
        match = regexp.search(self.ctx.stidy_output)
        self.ctx.space_group = int(match.group(1))

    def get_sites(self):
        regexp = re.compile(
            r'\s+([a-zA-Z]{1,2})(\d)+\s+([\w\(\)]{4,5})\s+([\d\/\.]+)\s+([\d\/\.]+)\s+([\d\/\.]+)\s+(\w+)\s+(\d+)')
        match_blocks = []
        for block in self.ctx.stidy_output.split('Wyckoff'):
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
        self.ctx.sites = match_blocks[0]

    def get_wyckoffs(self):
        wyckoffs = []
        for site in self.ctx.sites:
            wyckoffs.append(site[1].split('(')[-1].split(')')[0])
        self.ctx.wyckoffs = wyckoffs

    def get_cell(self):
        regexp = re.compile(r'^Cell :.*$', re.MULTILINE)
        match = regexp.search(self.ctx.stidy_output).group(0)
        abc = tuple(map(float, match.split(':')[-1].strip().split()[:3]))
        angles = tuple(map(float, match.split(':')[-1].strip().split()[3:]))
        self.ctx.cell = (abc, angles)

    def get_stidy_structure(self):
        cell = self.ctx.cell
        species = [site[-2] for site in self.ctx.sites]
        fractional_coords = [tuple(site[2:5]) for site in self.ctx.sites]
        lattice = Lattice.from_lengths_and_angles(abc=cell[0],
                                                  ang=cell[1])
        stidy_structure = Structure(lattice=lattice,
                                    species=species,
                                    coords=fractional_coords)
        self.ctx.stidy_structure = stidy_structure

    def get_fingerprint(self):
        ascii_iter = iter(ascii_uppercase)
        symbols = [site[5] for site in self.ctx.sites]
        symbol_map = {}
        for symbol in symbols:
            if symbol not in symbol_map.keys():
                symbol_map[symbol] = next(ascii_iter)
        wyckoffs = self.ctx.wyckoffs
        space_group = self.ctx.space_group
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
        self.out('fingerprint', Str(fingerprint))
        self.out('structure', StructureData(pymatgen=self.ctx.stidy_structure))
        print(len(self.ctx.structure.sites))
        print(len(self.ctx.stidy_structure.sites))
