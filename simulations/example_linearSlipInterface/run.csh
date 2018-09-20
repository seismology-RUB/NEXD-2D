#!/bin/tcsh
setenv OMP_NUM_THREADS $NSLOTS
echo "start program on"
echo $NSLOTS
mpirun -np $NSLOTS bin/solver
