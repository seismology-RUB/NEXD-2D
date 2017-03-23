#!/bin/tcsh
#qsub -l low -l arch=lx24-amd64 -R Y -cwd -hard -q "low.q@minos17" -l os=Debian_5\*_lenny -pe mpi-fu 1 run.csh
#qsub -l low -cwd -hard -q "low.q@minos18,low.q@minos19,low.q@minos20,low.q@minos21,low.q@minos22,low.q@minos23,low.q@minos24" -pe mpi-fu 96 run.csh
#qsub -l low -cwd -hard -q "low.q@minos20" -pe mpi-fu 16 run.csh
#qsub -l low -cwd -hard -q "low.q@minos26,low.q@minos27,low.q@minos23,low.q@minos24,low.q@minos25" -pe mpi-fu 128 run.csh
#qsub -l low -cwd -hard -q "low.q@minos26,low.q@minos27" -pe mpi-fu 64 run.csh
#qsub -l low -cwd -pe mpi-fu 128 run.csh
#qsub -l low -cwd -hard -q "low.q@minos26" -pe mpi-fu 2 run.csh
#setenv OMP_NUM_THREADS $NSLOTS
echo "start program on"
#echo $NSLOTS
mpirun -np 8 bin/solver
