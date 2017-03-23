#!/bin/bash
#

###################################################

# number of processes
NPROC=8

##################################################

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current directory
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p out

#rm -rf out/*

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



