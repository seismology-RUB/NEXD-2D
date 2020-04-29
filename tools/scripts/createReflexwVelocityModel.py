# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------
#   Copyright 2014-2020 Marc S. Boxberg (RWTH Aachen University, GER)
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
import sys
import numpy as np
import matplotlib.pyplot as plt

class mesh():

    def add_vp_values(self, filename):
        readvtk = open(filename, 'r')

        word = 'LOOKUP_TABLE'
        line = 'A A'
        while line.split()[0] != word:
            line = readvtk.readline()
            if line.split() == []:
                line = 'A A'
        for i in range(self.nelem):
            self.elem[i].add_vp(float(readvtk.readline()))
        self.type='vpvs'
        print('Added values for vp')

    def add_vs_values(self, filename):
        readvtk = open(filename, 'r')

        word = 'LOOKUP_TABLE'
        line = 'A A'
        while line.split()[0] != word:
            line = readvtk.readline()
            if line.split() == []:
                line = 'A A'
        for i in range(self.nelem):
            self.elem[i].add_vs(float(readvtk.readline()))
        self.type = 'vpvs'
        print('Added values for vs')

    def __init__(self, filename, type):
        self.nodes=[]
        self.elem=[]
        self.type=type

        if self.type in ['vp', 'vs']:
            pass
        else:
            sys.exit('STOP PROGRAM. Type has to be either vp or vs')

        readvtk=open(filename,'r')

        word = 'POINTS'
        line = 'A A'
        while line.split()[0] != word:
            line = readvtk.readline()
            if line.split()==[]:
                line='A A'
        self.nnodes = int(line.split()[1])
        for i in range(self.nnodes):
            coord = readvtk.readline().split()
            self.nodes.append(node(float(coord[0]), float(coord[2]), i))
        print('Initialized {:d} nodes'.format(self.nnodes))

        word = 'CELLS'
        line = 'A A'
        while line.split()[0] != word:
            line = readvtk.readline()
            if line.split()==[]:
                line='A A'
        self.nelem = int(line.split()[1])
        for i in range(self.nelem):
            elem_num = readvtk.readline().split()
            self.elem.append(elem(self.nodes[int(elem_num[1])], self.nodes[int(elem_num[2])], self.nodes[int(elem_num[3])], i))
        print('Initialized {:d} elements'.format(self.nelem))

        word = 'LOOKUP_TABLE'
        line = 'A A'
        while line.split()[0] != word:
            line = readvtk.readline()
            if line.split()==[]:
                line='A A'
        for i in range(self.nelem):
            self.elem[i].add_vp(float(readvtk.readline()))

        readvtk.close()

    def plot(self, type, show=True):
        x_plot = []
        y_plot = []
        vel_plot = []
        for elem in self.elem:
            x_plot.append(sum(elem.nodes[i].x for i in [0, 1, 2]) / 3.)
            y_plot.append(sum(elem.nodes[i].y for i in [0, 1, 2]) / 3.)
            if type == 'vp':
                vel_plot.append(elem.vp)
            else:
                vel_plot.append(elem.vs)
        plt.tricontourf(x_plot, y_plot, vel_plot)
        plt.colorbar()
        if show:
            plt.show()

    def make_profile(self, start, dir=[0,1], type='vp'):
        points_x=[]
        points_y=[]
        points_data=[]
        points=[]
        for elem in self.elem:
            intersec=[]
            for nodeset in [[0,1],[1,2],[2,0]]:
                # Test if parallel
                if not abs(dir[0]*(elem.nodes[nodeset[1]].y-elem.nodes[nodeset[0]].y)-dir[1]*(elem.nodes[nodeset[1]].x-elem.nodes[nodeset[0]].x)) < 1.E-5:
                    #Find paramter for intersection
                    b=((elem.nodes[nodeset[0]].x-start[0])*dir[1]-(elem.nodes[nodeset[0]].y-start[1])*dir[0])/((elem.nodes[nodeset[1]].y-elem.nodes[nodeset[0]].y)*dir[0]-(elem.nodes[nodeset[1]].x-elem.nodes[nodeset[0]].x)*dir[1])
                    if b > 0 and b < 1:
                        int_point=[elem.nodes[nodeset[0]].x+b*(elem.nodes[nodeset[1]].x-elem.nodes[nodeset[0]].x),elem.nodes[nodeset[0]].y+b*(elem.nodes[nodeset[1]].y-elem.nodes[nodeset[0]].y)]
                        intersec.append(int_point)
                        points.append(int_point)
            if len(intersec)==2:
                #points_x.append((intersec[0][0]+intersec[1][0])/2.)
                #points_y.append((intersec[0][1]+intersec[1][1])/2.)
                points_x.append(intersec[0][0])
                points_x.append(intersec[1][0])
                points_y.append(intersec[0][1])
                points_y.append(intersec[1][1])
                if type=='vp':
                    points_data.append(elem.vp)
                    points_data.append(elem.vp)
                else:
                    points_data.append(elem.vs)
                    points_data.append(elem.vs)
        return points_x, points_y, points_data, points


class node():
    def __init__(self, x, y, num):
        self.x = x
        self.y = y
        self.num = num

class elem():
    def __init__(self, node1, node2, node3, num):
        self.nodes = [node1, node2, node3]
        self.num = num
        self.vp = None
        self.vs = None

    def add_vp(self, vp):
        self.vp = vp

    def add_vs(self, vs):
        self.vs = vs


def write_vel_file(filename, position, depthlist, velocitylist):
    velfile = open(filename, 'w')
    velfile.write(' \n')
    velfile.write(' \n')
    velfile.write(' {0:9.4f} {1:9.4f}\n'.format(position[0], position[1]))
    for depth, velocity in zip(depthlist, velocitylist):
        velfile.write(' {0:9.4f} {1:9.4f}\n'.format(depth, velocity))
    velfile.close()


def sort_lists(list_depth, list_vel):
    depth_above = 0.0
    vel_above = 0.0

    depth, vel = zip(*sorted(zip(list_depth, list_vel)))
    depth = list(depth)
    vel = list(vel)
    
    for i in range(len(depth)):
        if i > 1 and depth[i] == depth_above and vel[i] != vel_above and vel_2above != vel_above:
            vel[i-1] = vel[i]
            vel[i] = vel_above
        depth_above = depth[i]
        vel_2above = vel_above
        vel_above = vel[i]

    return depth, vel
    

def reduce_lists(list_depth, list_vel):
    new_depth = []
    new_vel = []
    
    new_depth.append(list_depth[0])
    new_vel.append(list_vel[0])
    
    for i in range(1, len(list_depth)):
        if list_vel[i] != new_vel[-1]:
            new_depth.append(list_depth[i-1])
            new_vel.append(list_vel[i-1])
            new_depth.append(list_depth[i])
            new_vel.append(list_vel[i])
    
    if new_depth[-1] != list_depth[-1] or new_vel[-1] != list_vel[-1]:
        new_depth.append(list_depth[-1])
        new_vel.append(list_vel[-1])
    
    return new_depth, new_vel

def create_1d_velocity_model(grid, velfile, pos_x, type='vp', reducelist=True):
    print('create', end='', flush=True)
    _, list_pos_y, list_vel, _ = grid.make_profile([pos_x, 0.0], type=type)
    list_depth = abs(np.array(list_pos_y) - max(list_pos_y))
    list_depth, list_vel = sort_lists(list_depth, list_vel)
    if reducelist:
        list_depth, list_vel = reduce_lists(list_depth, list_vel)
    print(', write', end='', flush=True)
    write_vel_file(velfile, [pos_x, 0.0], list_depth, list_vel)
    print(', done.', flush=True)


def create_2d_model(outfiles):
    modelfile = open('reflexw_model.2DM', 'w')
    modelfile.write('{:4d}\n'.format(len(outfiles)))
    for outfile in outfiles:
        modelfile.write('{}\n'.format(outfile))
    modelfile.close()
    

if __name__ == '__main__':
    cl = argparse.ArgumentParser(description='Create 1D-velocity model (Reflexw) from 2D vtk velocity model. \n Please provide a list OR outfile and position.')
    cl.add_argument('infile', help="Specify the file with the velocity model (vtk).", type=str)
    cl.add_argument('--list_of_profiles', help="Provide a file with two columns, where the first contains the name(s) of the output files and the second column contains the location(s) of the profiles.", default=None)
    cl.add_argument('--outfile', help="Specify the location and name of the output file (vel).", type=str, default='out.VEL')
    cl.add_argument('--position', help="Specify the x-location of the 1D profile.", type=float, default=0.0)
    cl.add_argument('--type', help="Specify the velocity: vp (default) or vs",
                    type=str, choices=['vp', 'vs'], default='vp')
    args = cl.parse_args()

    print('Read vtk-file...')
    grid = mesh(args.infile, args.type)
    print('')
    if args.list_of_profiles is None:
        create_1d_velocity_model(grid, args.outfile, args.position, type=args.type)
    else:
        outfiles = []
        positions = []
        listfile = open(args.list_of_profiles, 'r')
        lines = listfile.readlines()
        for line in lines:
            outfile, position = line.split()
            if outfile[:-4].lower() != '.vel':
                outfile = outfile+'.VEL'
            outfiles.append(outfile)
            positions.append(float(position))
        for i, (out, pos) in enumerate(zip(outfiles, positions)):
            print('Profile {0:{2:d}d}/{1:d} at {3:9.4f} ({4}): '.format(i+1, len(outfiles), len(str(len(outfiles))), pos, out), end='', flush=True)
            create_1d_velocity_model(grid, out, pos, type=args.type)
        print('', flush=True)
        print('Write 2D Model...', end='', flush=True)
        create_2d_model(outfiles)
        print(' done.', flush=True)
