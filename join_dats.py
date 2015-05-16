#!/usr/bin/env python
# -*- coding: utf-8 -*-
"Join dat files"
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


def write_omegas(name, head, out):
    with open(name, 'w') as fout:
        if type(head[0]) == str:
            fout.write(head[0])
            for i in range(1, len(head)):
                fout.write("\t(%d)" % head[i])
            fout.write("\n")
        for i in out:
            fout.write('\t'.join(map(str, i)) + '\n')


ipat = argv[1]
rangt = eval(argv[2])
oupat = argv[3]
by = int(argv[4])
fcount = 0
for fname in range(*rangt):
    if fcount % by == 0:
        if fcount != 0:
            write_omegas(oupat % (fcount // by,), head, out)
        head = None
        out = []
    fcount += 1
    with open(ipat % fname) as fp:
        l = fp.readline()
        if l.startswith('#'):
            spl = l.split()
            lhead = [spl[0]]
            for i in spl[1:]:
                lhead.append(int(i[1:-1]))
            ctr = 0
        else:
            spl = l.split()
            lhead = [1] * len(spl)
            tline = list(map(float, spl))
            if head is not None:
                for i in range(1, len(tline)):
                    out[0][i] = (
                        out[0][i] * head[i] + tline[i]) / (head[i] + 1)
                ctr = 1
            else:
                out.append(tline)
        for l in fp:
            sp = l.split()
            if not sp:
                continue
            tline = []
            for f in sp:
                try:
                    tline.append(float(f))
                except ValueError:
                    tline.append(0)
            if head is not None:
                for i in range(1, len(tline)):
                    if lhead[i]:
                        out[ctr][i] = (out[ctr][i] * head[i] +
                                       tline[i] * lhead[i]) / \
                            (head[i] + lhead[i])
                ctr += 1
            else:
                out.append(tline)
        if head is not None:
            for i in range(1, len(head)):
                head[i] += lhead[i]
        else:
            head = lhead
ctr = fcount // by
# print fcount, by, fcount // by
if fcount % by > 0:
    ctr += 1
if head:
    write_omegas(oupat % (ctr,), head, out)
