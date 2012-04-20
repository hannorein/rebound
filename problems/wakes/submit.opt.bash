#!/bin/bash
#PBS -N wakes
#PBS -l nodes=22:ppn=12
#PBS -q t-A

cd ~/rebound/problems/wakes
mpiexec -np 264 ./nbody 3 >>log.txt
