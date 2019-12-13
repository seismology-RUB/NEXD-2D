#!/bin/bash
#

###################################################

# number of processes
NPROC=3

##################################################

echo "running example: `date`"
currentdir=`pwd`

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
mpirun -np $NPROC bin/solver
./bin/movie
