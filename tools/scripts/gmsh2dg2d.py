# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------
#   Copyright 2020 Marc S. Boxberg (RWTH Aachen University, GER)
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

import argparse
import meshio
import numpy as np
import os

class dg2d_mesh:
    def __init__(self, mesh: meshio.Mesh, folder='.'):
        print('Creating mesh object...', end=' ', flush=True)
        self.absorbfile = 'absorb'
        self.coordfile = 'coord'
        self.freefile = 'free'
        self.matfile = 'mat'
        self.matpropfile = 'matprop'
        self.meshfile = 'mesh'

        self.folder = folder
        self.mesh = mesh
        self.set_coord()
        self.set_mesh()
        self.set_matprop()
        self.set_mat()
        self.set_absorb()
        self.set_free()
        print('Done!')

    def set_coord(self):
        """
        Obtain the 2D coordinate matrix from the mesh.
        """
        self.coord = np.delete(self.mesh.points, 1, 1)

    def set_mesh(self):
        """
        Obtain the Cells as three coordinate IDs from the mesh.
        """
        if int(meshio.__version__.split('.')[0])>=4:
            self.cells = self.mesh.get_cells_type('triangle')+1
        else:
            self.cells = self.mesh.cells['triangle']+1

    def set_mat(self):
        """
        Obtain the link of material properties to the cells from the
        mesh.
        """
        if int(meshio.__version__.split('.')[0])>=4:
            self.mat = self.mesh.cell_data_dict['gmsh:physical']['triangle']
        else:
            self.mat = self.mesh.cell_data['triangle']['gmsh:physical']

    def set_matprop(self):
        """
        Obtain the material properties (physical groups) from the mesh.
        """
        self.matprop = dict()
        for physical_group in self.mesh.field_data.keys():
            properties = physical_group.split()
            if properties[0] in ['elastic']:
                self.matprop[self.mesh.field_data[physical_group][0]]={
                    'imat': 1,
                    'prop': str(1)+' '+list_to_string(properties[2:-1]),
                    'pml':  int(properties[-1]),
                    'flag': int(properties[1])}
            elif properties[0] in ['poro']:
                self.matprop[self.mesh.field_data[physical_group][0]]={
                    'imat': 2,
                    'prop': list_to_string(properties[2:]),
                    'pml':  0,
                    'flag': int(properties[1])}

    def set_absorb(self):
        """
        Obtain the nodes that belong to absorbing boundaries from the
        mesh.
        """
        self.absorb = []
        absorb_id = -1
        for (physical_group, ids) in self.mesh.field_data.items():
            if physical_group == 'absorb':
                absorb_id = ids[0]
        if absorb_id == -1: return
        if int(meshio.__version__.split('.')[0])>=4:
            celllist = self.mesh.cell_data_dict['gmsh:physical']['line']
        else:
            celllist = self.mesh.cell_data['line']['gmsh:physical']
        for i, cell in enumerate(celllist):
            if cell == absorb_id:
                if int(meshio.__version__.split('.')[0])>=4:
                    nodes = self.mesh.get_cells_type('line')[i]
                else:
                    nodes = self.mesh.cells['line'][i]
                for node in nodes:
                    if node not in self.absorb:
                        self.absorb.append(node)

    def set_free(self):
        """
        Obtain the nodes that belong to free surfaces / reflecting
        boundaries from the mesh.
        """
        self.free = []
        free_id = -1
        for (physical_group, ids) in self.mesh.field_data.items():
            if physical_group == 'free':
                free_id = ids[0]
        if free_id == -1: return
        if int(meshio.__version__.split('.')[0])>=4:
            celllist = self.mesh.cell_data_dict['gmsh:physical']['line']
        else:
            celllist = self.mesh.cell_data['line']['gmsh:physical']
        for i, cell in enumerate(celllist):
            if cell == free_id:
                if int(meshio.__version__.split('.')[0])>=4:
                    nodes = self.mesh.get_cells_type('line')[i]
                else:
                    nodes = self.mesh.cells['line'][i]
                for node in nodes:
                    if node not in self.free:
                        self.free.append(node)

    def write(self):
        """
        Method to write the files usable in NEXD 2D.
        """
        self._write_coord()
        self._write_mesh()
        self._write_matprop()
        self._write_mat()
        self._write_boundary(self.absorb, self.absorbfile)
        self._write_boundary(self.free, self.freefile)

    def _write_boundary(self, bdtype, bdname):
        """
        Writes the boundary file (either free or absorb) for NEXD 2D.
        """
        outfile = open(os.path.join(self.folder,bdname), 'w')
        outfile.write(str(len(bdtype))+'\n')
        i = -1
        for i, entry in enumerate(bdtype):
            print('Writing '+bdname+' (node '+str(i+1)+'/'+str(len(bdtype))+')...', end='\r', flush=True)
            outfile.write(str(entry+1)+'\n')
        outfile.close()
        print('Writing '+bdname+' (node '+str(i+1)+'/'+str(len(bdtype))+')... Done!')

    def _write_coord(self):
        """
        Writes the coordinate file for NEXD 2D.
        """
        outfile = open(os.path.join(self.folder,self.coordfile), 'w')
        outfile.write('%10i\n' % len(self.coord))
        for node, entry in enumerate(self.coord):
            print('Writing '+self.coordfile+' (node '+str(node+1)+'/'+str(len(self.coord))+')...', end='\r', flush=True)
            outfile.write(('%10i %20f %20f\n') % (node+1,entry[0],entry[1]))
        outfile.close()
        print('Writing '+self.coordfile+' (node '+str(node+1)+'/'+str(len(self.coord))+')... Done!')

    def _write_mat(self):
        """
        Writes the materials file for NEXD 2D.
        """
        outfile = open(os.path.join(self.folder,self.matfile), 'w')
        for i, mat in enumerate(self.mat):
            print('Writing '+self.matfile+' (element '+str(i+1)+'/'+str(len(self.mat))+')...', end='\r', flush=True)
            outfile.write(('%10i %10i %10i\n') % (i+1,self.matprop[mat]['flag'],self.matprop[mat]['pml']))
        outfile.close()
        print('Writing '+self.matfile+' (element '+str(i+1)+'/'+str(len(self.mat))+')... Done!')

    def _write_matprop(self):
        """
        Writes the material propierties file for NEXD 2D.
        """
        outfile = open(os.path.join(self.folder,self.matpropfile), 'w')
        outfile.write(str(len(self.matprop))+'\n')
        for i, (idnum, entry) in enumerate(self.matprop.items()):
            print('Writing '+self.matpropfile+' (material '+str(i+1)+'/'+str(len(self.matprop))+')...', end='\r', flush=True)
            outfile.write(str(i+1)+' '+entry['prop']+'\n')
        outfile.close()
        print('Writing '+self.matpropfile+' (material '+str(i+1)+'/'+str(len(self.matprop))+')... Done!')

    def _write_mesh(self):
        """
        Writes the mesh file for NEXD 2D.
        """
        outfile = open(os.path.join(self.folder,self.meshfile), 'w')
        outfile.write(str(len(self.cells))+'\n')
        for i, cell in enumerate(self.cells):
            print('Writing '+self.meshfile+' (element '+str(i+1)+'/'+str(len(self.cells))+')...', end='\r', flush=True)
            outfile.write(str(i+1)+' '+('%10i %10i %10i\n')% (cell[0],cell[1],cell[2]))
        outfile.close()
        print('Writing '+self.meshfile+' (element '+str(i+1)+'/'+str(len(self.cells))+')... Done!')

def list_to_string(list_):
    """
    Converts a list to a string with blanks as delimiters between the
    entries of the input list.

    :param list_: The list that has to be converted to a string.
    :type list_: list

    :returns: A string containing the elements of the input list.
    :rtype: str
    """
    final_string = str()
    for element in list_:
        final_string += str(element) + " "
    return final_string

def main(meshfile, folder='.'):
    """
    The main function that loads a given mesh and saves the files for
    NEXD 2D to a given folder.

    :param meshfile: String that gives the path to a (Gmsh-)file that
        contains the mesh that should be converted.
    :type meshfile: str
    :param folder: Output folder for the created NEXD 2D input files
        (default is '.').
    :type folder: str
    """
    print('Reading Gmsh file...', end=' ', flush=True)
    oldmesh = meshio.read(meshfile)
    print('Done!')
    newmesh = dg2d_mesh(oldmesh, folder)
    newmesh.write()
    print('Finished conversion to NEXD 2D.')

if __name__ == '__main__':
    cl = argparse.ArgumentParser(description = 'Convert gmsh-file (*.msh) to NEXD 2D format.')
    cl.add_argument('meshfile', help = "Specify the path to the meshfile to be converted.")
    cl.add_argument('--outputfolder', '-o', help = "Specify the output folder.", type=str, dest='folder', default='.')
    args = cl.parse_args()

    main(args.meshfile, args.folder)
