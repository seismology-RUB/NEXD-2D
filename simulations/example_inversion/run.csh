#!/bin/tcsh
echo "start program on"
echo $NSLOTS
mpirun -np $NSLOTS bin/solver -> 'out/output_solver'
