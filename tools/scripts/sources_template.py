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

#-------------------------------------------------------------
#   Module with templates for parfiles
#-------------------------------------------------------------
#  Template for gfdsv parfiles
#
def sources_head(nsrc):
    head_string = """# List of Source(s)
# Sources should be added sequentially. Otherwise problems will occur when reading the file

# Number of sources
nsrc        = {0:d}     # Number of sources

""".format(nsrc)
    return head_string

def sources_template(source, xsrc, zsrc, stf, f0, angle):
    """creates string containing parfile for gfdsv-computations"""
    parfile_string = """# Source {0:d} - Parameters
source      = {0:d}
xsource     = {1:f}     # x-value of the source
zsource     = {2:f}     # z-value of the source
delay       = 0.0       # time delay for the start of the source influence
sourcetype  = 0         # Type of source: 1 = Moment tensor, 0 = single force
stf         = {3:d}     # Type of source-time-function: 1 = green, 2 = ricker, 3 = sin^3, 4 = external
extwavelet  = data/1MHz0.5in_clean.txt #file with external wavelet; required if stf = 4
f0          = {4:f}     # center frequency of sft
factor      = 1e6       # factor of the sft
angle_force = {5:f}     # angle of the force action in the media
Mxx         = 0.0       # Momenttensor Mxx
Mzz         = 0.0       # Momenttensor Mzz
Mxz         = 1.0       # Momenttensor Mxz

""".format(source, xsrc, zsrc, stf, f0, angle)
    return parfile_string
