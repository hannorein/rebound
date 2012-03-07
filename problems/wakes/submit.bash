#!/bin/bash
#PBS -N wakes
#PBS -l nodes=6:ppn=48
#PBS -q t-A

cd ~/rebound/problems/wakes
mpiexec -np 288 ./nbody 10 >>log.txt
