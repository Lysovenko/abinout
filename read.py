#!/usr/bin/env python
# -*- coding: utf-8 -*-
"Read"
# wxRays (C) 2015 Serhii Lysovenko
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import os.path as osp
import numpy as np
from qfa import Cell


class Abinout_reader:
    def __init__(self, fname):
        if not osp.isfile(fname):
            raise IOError("file `%s' not found" % fname)
        self.the_fp = open(fname)
        self.out_num = -1

    def reset(self):
        self.the_fp.seek(0)
        self.out_num = -1

    def skip_to_next(self, delim='--- Iteration: '):
        line = self.the_fp.readline()
        while not line.startswith(delim):
            if line == '':
                raise IOError("no more OUTPUT items in the file")
            line = self.the_fp.readline()
        return line

    def read_floats(self, mask):
        sposition = self.the_fp.tell()
        line = self.the_fp.readline()
        while mask not in line:
            if line == '' or '--- Iteration: ' in line:
                self.the_fp.seek(sposition)
                return None
            line = self.the_fp.readline()
        proceed = True
        result = []
        while proceed:
            words = self.the_fp.readline().split()
            try:
                flts = tuple(map(float, words))
                result.append(flts)
            except ValueError:
                proceed = False
        self.the_fp.seek(sposition)
        return result

    def read_ccst(self):
        sposition = self.the_fp.tell()
        line = self.the_fp.readline()
        while not line.startswith(' Cartesian components of stress\
 tensor (hartree/bohr^3)'):
            if line == '' or '--OUTPUT--' in line:
                self.the_fp.seek(sposition)
                return None
            line = self.the_fp.readline()
        result = []
        for i in range(3):
            result.append(float(self.the_fp.readline()[14:30]))
        self.the_fp.seek(sposition)
        return result

    def read_mfloats(self, masks):
        sposition = self.the_fp.tell()
        results = {}
        line = self.the_fp.readline()
        while line.isspace():
            line = self.the_fp.readline()
        # print 'masks:', masks
        while not line.isspace():
            to_readline = True
            if line[2:4].isalpha():
                # print 'line:', line,
                for mask in masks:
                    if mask in line:
                        # print 'mask:', mask
                        sav_line = line
                        to_readline = False
                        proceed = True
                        result = []
                        while proceed:
                            line = self.the_fp.readline()
                            words = line.split()
                            try:
                                flts = np.array(map(float, words))
                                result.append(flts)
                            except ValueError:
                                proceed = False
                        if not result:
                            try:
                                result = float(sav_line.split()[-1])
                            except ValueError:
                                raise ValueError("mask '%s' is neither \
array or single float" % mask)
                        results[mask] = result
                        # print 'last line:', line
            if to_readline:
                line = self.the_fp.readline()
        return results

    def get_variables(self, vararr):
        self.the_fp.seek(0)
        self.skip_to_next(
            ' -outvars: echo values of preprocessed input variables --')
        curn = None
        curarr = None
        res = {}
        line = self.the_fp.readline()
        while not line.isspace():
            name = line[2:17].lstrip()
            if name != '':
                if curarr is not None:
                    res[curn] = curarr
                    curarr = None
                curn = name
                if curn in vararr:
                    curarr = []
            if curarr is not None:
                warr = line[17:].split()
                for i in warr:
                    try:
                        fv = float(i)
                    except ValueError:
                        continue
                    curarr.append(fv)
            line = self.the_fp.readline()
        if curarr is not None:
            res[curn] = curarr
        return res

    def __getitem__(self, item):
        if type(item) == str:
            return self.read_floats(item)
        elif type(item) == list:
            self.skip_to_next('---OUTPUT----------------------------------')
            return self.read_mfloats(item)
        else:
            try:
                self.skip_to_next()
            except IOError, err:
                raise IndexError(err)
            self.out_num += 1
            return self.out_num


class Qfa_cells:
    def __init__(self, reader):
        self.reader = reader
        var = reader.get_variables(['acell', 'amu', 'typat', 'natom'])
        typat = var['typat']
        self.acell = var['acell'][0] * 0.2645886054
        # 1 Bohr=0.5291772108 Angstroms
        typd = {}
        for i in map(int, typat):
            typd[i] = 1 + typd.get(i, 0)
        self.typat = []
        for i in typd:
            self.typat.append(typd[i])

    def __getitem__(self, item):
        self.reader[item]
        cell = Cell()
        cell.random(self.typat)
        cell.set_size(self.acell)
        try:
            for i, c in enumerate(self.reader['xred']):
                c = tuple([((x < 0. and x + 1. or x) - .5) * 2. for x in c])
                cell[i] = c
        except TypeError, err:
            raise IndexError(err)
        return cell


if __name__ == '__main__':
    from sys import stdout as sout
    sout.write('#hello world!\n')
    from sys import argv, stderr
    if len(argv) == 1:
        stderr.write('no input file\n')
        exit(1)
    reader = Abinout_reader(argv[1])
    var = reader.get_variables(['acell', 'amu', 'znucl', 'natom'])
    reader[0]
    vol = np.prod(var['acell'])
    for x in reader:
        P1 = reader.read_ccst()
        ress = reader[['vel', 'xcart', 'fcart', 'ekin', 'etot']]
        if ress == {}:
            print x
            continue
        vs = np.sum(ress['vel'], 0)
        P = np.dot(vs, vs)
        if P1 is None:
            P1 = 0
        else:
            P1 = - sum(P1) / 3. * 2.9421912e4
        Tx = ress['ekin'] * 2. * 105258.21342688377 / var['natom'][0]
        print '\t'.join(map(repr, [x, P, P1, Tx]))
    print '#', (var)
