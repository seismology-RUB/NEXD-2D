#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------
#   Copyright 2014-2018 Thomas Möller (Ruhr-Universität Bochum, GER)
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

#-----------------------------------------------------------------------
#   Creates the sources parameter file for the DG2D code from a template.
#   This is particular usefull when dealing with extended sources.
#   This script has to be executed in the location where the sources file has to be saved, or the needs to be copied to the desired location.
#----------------------------------------------------------------------
from sources_template import sources_template
from sources_template import sources_head
import argparse
import re
#-----------------------------------------------------------------------
#  set up commandline parsing
#
cl = argparse.ArgumentParser(description = 'Create sources file for DG2D code:')
cl.add_argument('src_file', help = "Name of sources parfile")
cl.add_argument('nsrc', help = "Total number of sources")
cl.add_argument('stf', help = "Select source-time-function: 1 = green, 2 = ricker, 3 = sin^3, 4 = external, 5 = sin^2")
cl.add_argument('xsrc', help = "X-coordinate of the first source")
cl.add_argument('zsrc', help = "Z-coordinate of the first source")
cl.add_argument('dx', help = "Increment between the sources along the x-axis")
cl.add_argument('dz', help = "Increment between the sources along the x-axis")
cl.add_argument('angle', help = "Angle of the force action.")
cl.add_argument('f0', help = "Central frequency of the stf")
args = cl.parse_args()

#Save parsed arguments
stf   = int(args.stf)
nsrc  = int(args.nsrc)
dx    = float(args.dx)
dz    = float(args.dz)
angle = float(args.angle)
f0 = float(args.f0)

#Create sources file
source_file = open(args.src_file, 'w')
#Generate Haeder
src_head = sources_head(nsrc)
source_file.write(src_head)

#write Sources
for source in range(nsrc):
    xsrc = float(args.xsrc) + (dx*source)
    zsrc = float(args.zsrc) + (dz*source)
    parfile_string = sources_template(source+1, xsrc, zsrc, stf, f0, angle)
    source_file.write(parfile_string)

