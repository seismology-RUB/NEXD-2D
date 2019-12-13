# -*- coding: utf-8 -*-
#!python
#!/usr/bin/env python

#-----------------------------------------------------------------------
#   Copyright 2014-2019 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
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

'''
This tool creates vtk files from any sorted nodeset in Trelis. This is
useful to show interfaces between different media in paraview while
looking at the wavefiels.
To create such vtk-file, simply create a nodeset in Trelis (1) and name
it "vtk TEXT" (2), where 'vtk' is the required keyword and 'TEXT' can be
anything upon your choice. 'TEXT' will also be used to name the vtk-file
as 'TEXT.vtk'.

Examples for commands in Trelis/Cubit:
(1) nodeset 3 node in curve 2
(2) nodeset 3 name "vtk interface_topSinemurian"

'''

import cubit
import vtk
import os.path as ospath
from sys import argv as sys_argv
from sys import exit as sys_exit


def main():
    try:
        inp_file = sys_argv[1]
        out_fold = sys_argv[2]
    except:
        print("USAGE: \n  createVtkOfInterfaceBetweenLayers.py input_file (.cub) output_folder \n e.g. \n createVtkOfInterfaceBetweenLayers.py cubit/testtri.cub out \n\n")
        sys_exit()

    start_cubit(inp_file)
    find_boundaries(out_fold)

def start_cubit(filename):
    cubit.init([""])
    cubit.cmd('set duplicate block elements on')
    command = 'open "{0}"'.format(filename)
    cubit.cmd(command)

def find_boundaries(outputfolder):
    vtk_found = False
    nodeset=cubit.get_nodeset_id_list()
    for n in nodeset:
        name=cubit.get_exodus_entity_name('nodeset',n)
        if name.find('vtk') >=0:
            vtk_found = True
            print('found interface nodes to create vtk')
            values=name.split(" ")
            title = values[1]
            nodes_vtk=cubit.get_nodeset_nodes(n)
            write_vtk_file(nodes_vtk,title,outputfolder)
    if not vtk_found:
        print('no interface nodes to create vtk found')


def write_vtk_file(nodeset,title,folder):
    writer = vtk.vtkUnstructuredGridWriter()
    filename = ospath.join(folder,'{0}.vtk'.format(title))
    writer.SetFileName(filename)
    writer.SetInputData(MakePolyLine(nodeset))
    print('Write {0}'.format(filename))
    writer.Write()


def MakePolyLine(nodeset):
    # A polyline is a cell that represents a set of 1D lines
    numberOfVertices = len(nodeset)

    points = vtk.vtkPoints()
    for node in nodeset:
        x,y,z=cubit.get_nodal_coordinates(node)
        points.InsertNextPoint(x, z, y)

    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(numberOfVertices)

    for i in range(0, numberOfVertices):
        polyline.GetPointIds().SetId(i, i)

    ug = vtk.vtkUnstructuredGrid()
    ug.SetPoints(points)
    ug.InsertNextCell(polyline.GetCellType(), polyline.GetPointIds())

    return ug


if __name__ == '__main__':
    main()
