#!/usr/bin/env python
# -*- coding: utf-8 -*-
"Diffusion"
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

from sys import argv
from read import Abinout_reader, Qfa_cells

in_name = argv[1]
skip = int(argv[2])
reader = Abinout_reader(in_name)
dt = reader.get_variables(('dtion')).get('dtion', (100.,))[0] * 2.418884e-5
# picoseconds
reader.reset()
cells = Qfa_cells(reader)
for i in xrange(skip):
    refc = cells[i]
types = len(refc.atoms_of_types())
for i, cell in enumerate(cells):
    sh = [cell.mid_r2_shift(refc, t) for t in range(types)]
    print("%g\t%s" % (i*dt, '\t'.join(map(str, sh))))
