#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
from obspy import Trace, UTCDateTime

def conv2sac(root, component, datatype, conv_format):

    fpattern = 'seismo.' + component + '.*.' + 'sd' + datatype

    chan_conv = dict(X='1', Y='2', Z='3', R='R', T='T')

    files = glob.glob1(root, fpattern)

    for fn in files:
        data = np.loadtxt(fn)
        chan, station = fn.split('.')[1].strip().upper(), fn.split('.')[2].strip()

        time = data[:,0]
        stime = UTCDateTime(time[0])
        delta = time[1] - time[0]
        freq = 1 / delta
        data = data[:,1]
        tr = Trace()
        tr.data = data
        tr.stats.station = station
        tr.stats.channel = chan_conv[chan]
        tr.stats.starttime = stime
        tr.stats.delta = delta
        tr.stats.sampling_rate = freq

        outfile = '{0}.{1}.{2}'.format(station, tr.stats.channel, conv_format.lower())
        tr.write(outfile, format = conv_format)

if __name__ == '__main__':
    cl = argparse.ArgumentParser(description = 'Create sac-files from seismogram data')
    cl.add_argument('dir', help = "Specify the directory of the files to be converted")
    cl.add_argument('comp', help = "Specify the component of the seismogram", type=str, choices = ['x', 'y', 'z', 'r', 't'])
    cl.add_argument('type', help = "Specify the seismogram: Velocity (v), Acceleration (a), Displacement (u)")
    cl.add_argument('--format', '-f', help = "Specify the output format for the written seismogram files.", type=str, dest='format', default='SAC')
    args = cl.parse_args()

    conv2sac(args.dir, args.comp, args.type, args.format)
