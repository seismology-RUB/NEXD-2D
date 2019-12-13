# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------
#   Copyright 2019 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
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

class Source:
    def __init__(self, dgfile, number):
        self.dgfile = dgfile
        
        # global attributes
        self.source = number
        self._f0 = 0.
        self._factor = 0.
        self._Mxx = 0.
        self._Mzz = 0.
        self._Mxz = 0.
        
        # SPECFEM2D attributes
        self._source_surf = False
        self._xs = 0.0
        self._zs = 0.0
        self._source_type = 0
        self._time_function_type = 0
        self._name_of_source_file = ''
        self._burst_band_width = 0.
        self._tshift = 0.
        self._anglesource = 0.
        
        # DG2D attributes
        self._xsource = 0.0
        self._zsource = 0.0
        self._delay = 0.0
        self._sourcetype = 0
        self._stf = 0
        self._extwavelet = ''
        self._angle_force = 0.0
    
    def set_f0(self, f0):
        self._f0 = float(f0)
    
    def set_factor(self, factor):
        self._factor = float(factor)
    
    def set_Mxx(self, Mxx):
        self._Mxx = float(Mxx)
        
    def set_Mzz(self, Mzz):
        self._Mzz = float(Mzz)
    
    def set_Mxz(self, Mxz):
        self._Mxz = float(Mxz)
    
    def set_xs(self, xs):
        self._xs = float(xs)
        self._xsource = self._xs
        
    def set_zs(self, zs):
        self._zs = float(zs)
        self._zsource = self._zs
    
    def set_source_type(self, source_type):
        self._source_type = int(source_type)
        self._sourcetype = int(source_type) - 1
        
    def set_time_function_type(self, time_function_type):
        translation = {1:2,     # Ricker
                       2:None,  # first derivative
                       3:1,     # Gauss
                       4:None,  # Dirac
                       5:None,  # Heaviside
                       8:4,     # External
                       9:None   # Burst
                      }
        if translation[int(time_function_type)] is None:
            raise NotImplementedError('time_function_type {:d} not available in dg2d.'.format(time_function_type))
        else:
            self._time_function_type = int(time_function_type)
            self._stf = translation[self._time_function_type]
    
    def set_name_of_source_file(self, name_of_source_file):
        self._name_of_source_file = name_of_source_file
        self._extwavelet = name_of_source_file
    
    def set_burst_band_width(self, burst_band_width):
        self._burst_band_width = float(burst_band_width)
        # NOT IMPLEMENTED IN DG2D
    
    def set_source_surf(self, source_surf):
        translate = {'.true.': True,
                     '.false.': False}
        self._source_surf = translate[source_surf]
        # NOT NECESSARY IN DG2D
    
    def set_tshift(self, tshift):
        self._tshift = float(tshift)
        self._delay = self._tshift
    
    def set_anglesource(self, anglesource):
        self._anglesource = float(anglesource)
        self._angle_force = self._anglesource
    
    def write(self):
        self.dgfile.write('# Source {} - Parameters\n'.format(self.source))
        self.dgfile.write('source      = {}  # Running Number of the sources\n'.format(self.source))
        self.dgfile.write('xsource     = {}  # x-value of the source\n'.format(self._xsource))
        self.dgfile.write('zsource     = {}  # z-value of the source\n'.format(self._zsource))
        self.dgfile.write('delay       = {}  # time delay for the source\n'.format(self._delay))
        self.dgfile.write('sourcetype  = {}  # Type of source: 1 = Moment tensor, 0 = single force\n'.format(self._sourcetype))
        self.dgfile.write('stf         = {}  # Type of source-time-function: 1 = gauss, 2 = ricker, 3 = sin^3, 4 = external\n'.format(self._stf))
        self.dgfile.write('extwavelet  = {}  #file with external wavelet; required if stf = 4\n'.format(self._extwavelet))
        self.dgfile.write('f0          = {}  # center frequency of stf\n'.format(self._f0))
        self.dgfile.write('factor      = {}  # factor of the stf\n'.format(self._factor))
        self.dgfile.write('angle_force = {}  # angle of the force action in the media\n'.format(self._angle_force))
        self.dgfile.write('Mxx         = {}  # Momenttensor Mxx\n'.format(self._Mxx))
        self.dgfile.write('Mzz         = {}  # Momenttensor Mzz\n'.format(self._Mzz))
        self.dgfile.write('Mxz         = {}  # Momenttensor Mxz\n'.format(self._Mxz))
        self.dgfile.write('\n')

class Statistics:
    def __init__(self):
        self._f0 = {'min': None, 'max':None}
        self._xsource = {'min': None, 'max':None}
        self._zsource = {'min': None, 'max':None}
        self._delay = {'min': None, 'max':None}
    
    def check(self, src):
        self.set_f0(src._f0)
        self.set_delay(src._delay)
        self.set_xsource(src._xsource)
        self.set_zsource(src._zsource)
    
    def set_f0(self, value):
        if self._f0['min']:
            self._f0['min'] = min(self._f0['min'], value)
        else:
            self._f0['min'] = value
        if self._f0['max']:
            self._f0['max'] = max(self._f0['max'], value)
        else:
            self._f0['max'] = value
        
    def set_xsource(self, value):
        if self._xsource['min']:
            self._xsource['min'] = min(self._xsource['min'], value)
        else:
            self._xsource['min'] = value
        if self._xsource['max']:
            self._xsource['max'] = max(self._xsource['max'], value)
        else:
            self._xsource['max'] = value
        
    def set_zsource(self, value):
        if self._zsource['min']:
            self._zsource['min'] = min(self._zsource['min'], value)
        else:
            self._zsource['min'] = value
        if self._zsource['max']:
            self._zsource['max'] = max(self._zsource['max'], value)
        else:
            self._zsource['max'] = value
        
    def set_delay(self, value):
        if self._delay['min']:
            self._delay['min'] = min(self._delay['min'], value)
        else:
            self._delay['min'] = value
        if self._delay['max']:
            self._delay['max'] = max(self._delay['max'], value)
        else:
            self._delay['max'] = value
    
    def __str__(self):
        return 'xsource = [{},{}], zsource = [{},{}], f0 = [{},{}], delay = [{},{}]'.format(
                self._xsource['min'],self._xsource['max'],
                self._zsource['min'],self._zsource['max'],
                self._f0['min'],self._f0['max'],
                self._delay['min'],self._delay['max'])

def start_new_source(source, dgfile, number, stats):
    if source:
        stats.check(source)
        source.write()
    source = Source(dgfile, number)
    return source

def convert(par_file, filename_specfem2d, filename_dg2d):
    source = None
    stats = Statistics()
    nsources = 0
    
    specfemparfile = open(par_file, 'r')
    line = specfemparfile.readline()
    while line:
        content = line.split()
        if content and content[0] == 'NSOURCES':
            nsources = int(content[2])
            break
        line = specfemparfile.readline()
    if nsources == 0:
        raise ValueError('Parameter NSOURCES not found in Par_file {}'.format(specfemparfile))
    
    dgfile = open(filename_dg2d, 'w')
    dgfile.write('# List of Source(s)\n')
    dgfile.write('# Sources should be added sequentially. Otherwise problems will occur when reading the file\n')
    dgfile.write('\n')
    dgfile.write('# Number of sources\n')
    dgfile.write('nsrc        = {:d}       # Total number of sources\n'.format(nsources))
    dgfile.write('\n')
    
    specfemfile = open(filename_specfem2d, 'r')
    line = specfemfile.readline()

    while line:
        content = line.split()
        if content[0] == '#source':
            source = start_new_source(source, dgfile, int(content[1]), stats)
        else:
            handler = {'source_surf' : source.set_source_surf,
               'xs' : source.set_xs,
               'zs' : source.set_zs,
               'source_type' : source.set_source_type,
               'time_function_type' : source.set_time_function_type,
               'name_of_source_file' : source.set_name_of_source_file,
               'burst_band_width' : source.set_burst_band_width, 
               'f0' : source.set_f0,
               'tshift' : source.set_tshift,
               'anglesource' : source.set_anglesource,
               'Mxx' : source.set_Mxx,
               'Mzz' : source.set_Mzz,
               'Mxz' : source.set_Mxz,
               'factor' : source.set_factor
               }
            handler[content[0]](content[2])
        line = specfemfile.readline()

    stats.check(source)
    print(stats)
    
    source.write()
    
    specfemfile.close()
    dgfile.close()

if __name__ == '__main__':
    cl = argparse.ArgumentParser(description='Convert SPECFEM2D source file to DG2D source file.')
    cl.add_argument('specfem2d_Par_file', help="Specify the SPECFEM2D Par_file", type=str)
    cl.add_argument('specfem2d_source_file', help="Specify the SPECFEM2D source file", type=str)
    cl.add_argument('dg2d_source_file', help="Specify the DG2D source file", type=str)
    args = cl.parse_args()

    convert(args.specfem2d_Par_file, args.specfem2d_source_file, args.dg2d_source_file)
