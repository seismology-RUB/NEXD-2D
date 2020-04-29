# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------
#   Copyright 2015-2020 Andre Lamert (Ruhr-Universität Bochum, GER)
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
# -----------------------------------------------------------------------

import os
import shutil

pathin='out/'
pathout='inversion/'

for files in os.listdir(pathin):
    if files.endswith('.sdu'):
        slash=files.rfind('/')
        newname=pathout+files
        shutil.copy(pathin+files,pathout+files+'.002.meas')
        print files, '   ', newname
