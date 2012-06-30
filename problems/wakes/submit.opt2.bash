#!/bin/bash
#PBS -o log.out
#PBS -e log.err
#PBS -N wakes
#PBS -l nodes=4:ppn=48
#PBS -q t-A
module purge
module load mvapich2/gcc

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=48
export MV2_ENABLE_AFFINITY=0

# Number 1
#mpiexec -npernode 1 ./nbody     
# Number 2
mpiexec -npernode 1 ./nbody --sigma=2400 --length=1000
# Number 3
#mpiexec -npernode 1 ./nbody --rmin=0.256 --rmax=2.56 

# Number 2 with twice the optical depth
#mpiexec -npernode 1 ./nbody --sigma=4800 --length=2000

# Number 2 with twice the optical depth and a size distribution
#mpiexec -npernode 1 ./nbody --sigma=4800 --length=2000 --rmin=0.5

