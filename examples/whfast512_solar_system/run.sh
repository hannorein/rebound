#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --time=1:00:00
#SBATCH --job-name=rebound_benchmark
#SBATCH --output=/scratch/r/rein/rein/output.txt
#SBATCH --mail-type=FAIL
 
rsync -ua --progress --exclude=".*" $HOME/git/rebound $SCRATCH/rebound 

cd $SCRATCH/rebound/rebound/examples/avx512_performance

module load intel
make clean
make -j 40

rm output100_*.txt

 
for gr in {0..0}
do
    for faster in {0..1}
    do
        for N in {1..8}
        do
            for i in {0..79}
            do
                ./rebound $gr $faster $N $i &
            done
            wait
        done
    done
done

module unload intel
module load gcc
make clean
make -j 40

 
for gr in {0..0}
do
    for faster in {0..1}
    do
        for N in {1..8}
        do
            for i in {0..79}
            do
                ./rebound $gr $faster $N $i &
            done
            wait
        done
    done
done

