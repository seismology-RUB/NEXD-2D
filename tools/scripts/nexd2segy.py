# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------
#   Copyright 2019 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
#   Copyright 2019 Kaan Cökerim (Ruhr-Universität Bochum, GER)
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
import os
import sys

import numpy as np
from obspy import Trace, Stream, UTCDateTime
from obspy.core import AttribDict
from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader

trace_identification_code = dict(x=14, z=12, p=11)
trace_value_measurement_unit = dict(u=5, v=6, a=7, p=1)


class Source:
    def __init__(self, file):
        self._angle_force = 0.
        self._x_source = 0.
        self._z_source = 0.
        self._source_type_orientation = -1

        self._read(file)

    def __str__(self):
        return 'Source at x={:.3f}, z={:.3f} with an angle of {:.3f}.'.format(self._x_source, self._z_source,
                                                                              self._angle_force)

    def _read(self, file):
        keywords = {'xsource': self.__set_x_source,
                    'zsource': self.__set_z_source,
                    'angle_force': self.__set_angle_force}
        with open(file, 'r') as source_file:
            for line in source_file:
                key = line.split(' ')[0]
                if key in keywords:
                    keywords[key](line.split()[2])

    @property
    def angle_source(self):
        return self._angle_force

    @angle_source.setter
    def angle_source(self, value):
        self.__set_angle_force(value)

    @property
    def source_type_orientation(self):
        return self._source_type_orientation

    @property
    def x_source(self):
        return self._x_source

    @x_source.setter
    def x_source(self, value):
        self.__set_x_source(value)

    @property
    def z_source(self):
        return self._z_source

    @z_source.setter
    def z_source(self, value):
        self.__set_z_source(value)

    def __set_angle_force(self, value):
        self._angle_force = float(value)
        if self._angle_force in [0.0, 180.0]:
            self._source_type_orientation = 4
        elif self._angle_force in [90.0, 270.0]:
            self._source_type_orientation = 6

    def __set_x_source(self, value):
        self._x_source = float(value)

    def __set_z_source(self, value):
        self._z_source = float(value)


def read_stations_file(file):
    with open(file, 'r') as stations_file:
        stations_list = []
        number_of_stations = 0
        for i, line in enumerate(stations_file):
            if i == 0:
                number_of_stations = int(line.split(' ')[2])
            else:
                try:
                    stations_list.append([float(line.split()[i]) for i in range(1, 3)])
                except:
                    continue
        if not (len(stations_list) == number_of_stations):
            raise ValueError('nrec not equal to number of given station coordinates.')
        return stations_list


def read_nexd_2d_seismogram(root, source, x_rec, z_rec, component, data_type, number, shotpoint_number, scaling=-100,
                            starttime=UTCDateTime(2019, 1, 1, 0, 0, 0, 0), resampling_rate=None, SNR=None):
    if scaling > 0:
        scale = 1 / scaling
    elif scaling < 0:
        scale = -scaling
    else:
        raise ValueError('scaling must not be equal to zero.')

    file = os.path.join(root, 'seismo.{0}.{1:07d}.sd{2}'.format(component, number, data_type))
    file_data = np.genfromtxt(file, unpack=True)
    tr = Trace()
    tr.data = file_data[1]

    tr.stats.starttime = starttime
    tr.stats.delta = file_data[0][1] - file_data[0][0]

    if resampling_rate:
        tr.resample(resampling_rate)

    if tr.stats.npts > 32767:
        factor = tr.stats.npts / 32767.
        print('resample (too many points) to {:f}'.format(tr.stats.sampling_rate / factor))
        tr.resample(tr.stats.sampling_rate / factor)

    if SNR:
        #creating array with signal energies
        Psig = tr.data**2
        
        # calculating average signal power as the root mean square of the energy array
        avg_Psig = np.mean(Psig)
        
        #  Pnoise = Psig / SNR
        avg_Pnoise = avg_Psig / SNR
        # creating normal distributed noise sample:
        # For a Gaussian random variable X, the average power E[X^2] is E[X^2] = mu^2 + sigma^2
        # mu = mean; sigma = standard deviation
        # for white noise, mu = 0 and the average power is then equal to the variance sigma^2.
        mean_noise = 0
        noise = np.random.normal(mean_noise, np.sqrt(avg_Pnoise), tr.stats.npts) # 1: mean, 2: standard deviation, 3: length of noise array

        # adding noise to pure signal
        tr.data += noise

    tr.data = np.require(tr.data, dtype=np.float32)

    if not hasattr(tr.stats, 'segy.trace_header'):
        tr.stats.segy = {}
        tr.stats.segy.trace_header = SEGYTraceHeader()
    tr.stats.segy.trace_header.coordinate_units = 1
    tr.stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group = int(
        np.abs(x_rec - source.x_source) * scale)
    tr.stats.segy.trace_header.group_coordinate_x = int(x_rec * scale)
    tr.stats.segy.trace_header.receiver_group_elevation = int(z_rec * scale)
    tr.stats.segy.trace_header.shotpoint_number = shotpoint_number
    tr.stats.segy.trace_header.scalar_to_be_applied_to_all_coordinates = scaling
    tr.stats.segy.trace_header.scalar_to_be_applied_to_all_elevations_and_depths = scaling
    tr.stats.segy.trace_header.source_coordinate_x = int(source.x_source * scale)
    tr.stats.segy.trace_header.source_type_orientation = source.source_type_orientation
    tr.stats.segy.trace_header.surface_elevation_at_source = int(source.z_source * scale)
    tr.stats.segy.trace_header.trace_identification_code = trace_identification_code[component]
    tr.stats.segy.trace_header.trace_sequence_number_within_line = number
    tr.stats.segy.trace_header.trace_value_measurement_unit = trace_value_measurement_unit[data_type]
    tr.stats.segy.trace_header.x_coordinate_of_ensemble_position_of_this_trace = int(x_rec * scale)

    return tr


def read_nexd_2d_source(root, source, scaling=-100, starttime=UTCDateTime(2019, 1, 1, 0, 0, 0, 0), resampling_rate=None):
    if scaling > 0:
        scale = 1 / scaling
    elif scaling < 0:
        scale = -scaling
    else:
        raise ValueError('scaling must not be equal to zero.')
    
    file_data = np.genfromtxt(root, unpack=True)
    tr = Trace()
    tr.data = file_data[1]

    tr.stats.starttime = starttime
    tr.stats.delta = file_data[0][1] - file_data[0][0]

    if resampling_rate:
        tr.resample(resampling_rate)

    if tr.stats.npts > 32767:
        factor = tr.stats.npts / 32767.
        tr.resample(tr.stats.sampling_rate / factor)

    tr.data = np.require(tr.data, dtype=np.float32)
    
    if not hasattr(tr.stats, 'segy.trace_header'):
        tr.stats.segy = {}
        tr.stats.segy.trace_header = SEGYTraceHeader()
    tr.stats.segy.trace_header.coordinate_units = 1
    tr.stats.segy.trace_header.scalar_to_be_applied_to_all_coordinates = scaling
    tr.stats.segy.trace_header.scalar_to_be_applied_to_all_elevations_and_depths = scaling
    tr.stats.segy.trace_header.source_coordinate_x = int(source.x_source * scale)
    tr.stats.segy.trace_header.source_type_orientation = source.source_type_orientation
    tr.stats.segy.trace_header.surface_elevation_at_source = int(source.z_source * scale)
    tr.stats.segy.trace_header.trace_identification_code = 6
    tr.stats.segy.trace_header.trace_value_measurement_unit = trace_value_measurement_unit['v']

    return tr

def read_nexd_2d_survey(root, out, component, data_type, shotpoint_number=1, resampling_rate=None, sourcefile=None, SNR=None):
    print('read source file...')
    source = Source(file=os.path.join(root, 'data', 'source'))
    print(source)
    print('read stations file...')
    stations = read_stations_file(file=os.path.join(root, 'data', 'stations'))

    st = Stream()
    if sourcefile:
        print('read sourcefile...')
        st.append(read_nexd_2d_source(root=sourcefile, source=source, resampling_rate=resampling_rate))
    
    for i, station in enumerate(stations):
        station_id = i + 1
        print('read station {0:{4:d}d}/{1:d} at x={2:.3f}, z={3:.3f}...'.format(station_id, len(stations), station[0],
                                                                                station[1], len(str(len(stations)))))
        st.append(read_nexd_2d_seismogram(root=os.path.join(root, 'out'), source=source, x_rec=station[0],
                                          z_rec=station[1], component=component, data_type=data_type, number=station_id,
                                          shotpoint_number=shotpoint_number, resampling_rate=resampling_rate, SNR=SNR))

    st.stats = AttribDict()
    st.stats.textual_file_header = 'Textual Header!'
    st.stats.binary_file_header = SEGYBinaryFileHeader()
    st.stats.binary_file_header.trace_sorting_code = 5  # Common source point

    print('write segy file...')
    st.write(out, format='SEGY', data_encoding=1, byteorder=sys.byteorder)


if __name__ == '__main__':
    cl = argparse.ArgumentParser(description='Create SEGY-files from seismograms')
    cl.add_argument('root_dir', help="Specify the main directory of the simulation.")
    cl.add_argument('output_file', help="Specify the location and name of the output file.")
    cl.add_argument('--component', help="Specify the component of the seismogram (default = 'z').", type=str,
                    choices=['x', 'y', 'z', 'p', 'r', 't'], default='z')
    cl.add_argument('--type',
                    help="Specify the seismogram: Velocity (v), Acceleration (a), Displacement (u) (default = 'v').",
                    type=str, choices=['v', 'a', 'u'], default='v')
    cl.add_argument('--shotpoint_number', help="Give the number of the shotpoint.", type=int, default=1)
    cl.add_argument('--resampling_rate', help="Give the sampling rate for resampling. Default is None.", type=float, default=None)
    cl.add_argument('--reference_trace', help="Provide the source wavelet (at the same sampling rate as the seismograms!).",
                    type=str, default=None)
    cl.add_argument('--SNR', help="Add artificial white noise with this signal-to-noise ratio. Default is None.", type=float, default=None)
    args = cl.parse_args()

    read_nexd_2d_survey(args.root_dir, args.output_file, component=args.component, data_type=args.type,
                        shotpoint_number=args.shotpoint_number, resampling_rate=args.resampling_rate,
                        sourcefile=args.reference_trace, SNR=args.SNR)

