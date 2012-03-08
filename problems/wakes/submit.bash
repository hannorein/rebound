#!/bin/bash
#PBS -N wakes
#PBS -l nodes=12:ppn=48
#PBS -q t-A

cd ~/rebound/problems/wakes
mpiexec -np 576 ./nbody 5 >>log.txt
