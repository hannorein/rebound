#!/bin/bash
#PBS -N wakes
#PBS -l nodes=1:ppn=6,walltime=24:00:00
#PBS -m n

module load mpt/2.05
cd $PBS_O_WORKDIR
export NPROCS=`wc -l <$PBS_NODEFILE`

mpirun -np $NPROCS ./nbody --boxsize=10 --root_nx=6 --root_ny=6
