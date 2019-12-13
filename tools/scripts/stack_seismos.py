#!/usr/bin/env python

#-----------------------------------------------------------------------
#   Copyright 2014-2019 Thomas Möller (Ruhr-Universität Bochum, GER)
#
#   This file is part of NEXD 2D.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with NEXD 2D. If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------

import argparse
import numpy as np
import glob

def stack_seismos(root, component, datatype):
    fpattern = 'seismo.' + component + '.*.' + 'sd' + datatype

    chan_conv = dict(X='1', Y='2', Z='3', R='R', T='T')

    files = glob.glob1(root, fpattern)

    outfile = 'seismo.' + component + '.stacked.sds' + datatype

    output = open(outfile, 'w', 0)
    i = 0
    for file in files:
        data = np.loadtxt(file)
        if i == 0:
            stack = data[:,1]
            time = data[:,0]
        else:
            new_stack = data[:,1]
            stack += new_stack
        i += 1
    stack = stack/i

    for item in range(data.shape[0]):
        output.write(str(time[item]) + "   " + str(stack[item]) + "\n")

    output.close()

if __name__ == '__main__':
    cl = argparse.ArgumentParser(description = 'Create a stacked file from given seismograms')
    cl.add_argument('dir', help = "Specify the directory of the files to be converted")
    cl.add_argument('comp', help = "Specify the component of the seismogram", type=str, choices = ['x', 'y', 'z', 'r', 't'])
    cl.add_argument('type', help = "Specify the seismogram: Velocity (v), Acceleration (a), Displacement (u)")
    args = cl.parse_args()

    stack_seismos(args.dir, args.comp, args.type)
