#!/bin/bash
#
echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current directory
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p out
mkdir -p adjoint
mkdir -p inversion

rm -f run.csh.*
rm -rf out/*
rm -rf cubit/matpropiter*
rm -rf adjoint/*
rm -rf inversion/seismo*src*
rm -rf inversion/seismo*filter*
rm -rf inversion/*bin
rm -rf inversion/*vtk
rm -rf inversion/*temp*
rm -rf inversion/*run*

# compiles executables in root directory
cd ../../
#make clean
make -j 4 all > tmp.log
cd $currentdir

# copy program
cd bin/
cp ../../../bin/mesher .
cp ../../../bin/solver .
cp ../../../bin/movie .
cd ../

# mesh
./bin/mesher

# stores setup
cp data/parfile out/
cp data/stations out/
cp data/source out/



