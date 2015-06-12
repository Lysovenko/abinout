#!/usr/bin/env python
# -*- coding: utf-8 -*-
"Extract cells"
# abinout (C) 2015 Serhii Lysovenko
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

from read import Abinout_reader
from sys import argv, stderr
from qfa import Cell, RDD
in_name = argv[1]
pref = argv[2]


def xred2c(x):
    if x < 0.:
        x += 1.
    if x > 1.:
        x -= 1.
    return (x - .5) * 2.

reader = Abinout_reader(in_name)
var = reader.get_variables(["acell", "amu", "typat", "natom"])
typat = var["typat"]
typd = {}
for i in map(int, typat):
    typd[i] = 1 + typd.get(i, 0)
typat = []
for i in range(1, len(typd) + 1):
    typat.append(typd[i])
print typat
reader[0]
for x in reader:
    cel = Cell()
    cel.random(typat)
    cel.set_size(var["acell"][0] * 0.2645886096)
    for i, c in enumerate(reader["xred"]):
        c = tuple(map(xred2c, c))
        cel[i] = c
    if '-rdf' in argv:
        rdd = RDD(cel, .05, .05, 0., 0)
        rdd.calculate_rdf()
        rdd.save(pref % x)
    else:
        cel.save(pref % x)
