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

rm -rf out/*

# compiles executables in root directory
echo "Compiling program..."
cd ../../
#make clean
make -j 4 all #> tmp.log
if [ $? -ne 0 ]
then
    echo "make failed. Abort..."
    cd $currentdir
    exit 1
fi
cd $currentdir

# copy program
echo "Copying binaries..."
cd bin/
cp ../../../bin/mesher .
cp ../../../bin/solver .
cp ../../../bin/movie .
cd ../

# stores setup
echo "Storing input parameters..."
cp data/parfile out/
cp data/stations out/
cp data/source out/

# mesh
echo "Execute mesher..."
./bin/mesher

echo "done at: `date`"