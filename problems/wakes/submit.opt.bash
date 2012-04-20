#!/bin/bash
#PBS -N wakes
#PBS -l nodes=1:ppn=6,walltime=48:00:00
#PBS -q n-A

MPI=/usr/local/mvapich2/mvapich2.1.8rc1-gcc
PATH=${MPI}/bin:${PATH}
LD_LIBRARY_PATH=${MPI}/lib

cd ~/rebound/problems/wakes
mpirun ./nbody --boxsize=10 --root_nx=6 --root_ny=6

