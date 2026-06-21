#!/bin/bash
#SBATCH --account=def-
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=6:00:00
#SBATCH --job-name=wf512energy
#SBATCH --output=ensemble_%j.log

set -u
export EPREFIX=/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3
export PATH="$EPREFIX/usr/bin:$PATH"
cd /lustre09/project/6098391/rishit/rebound/examples/whfast512_energy

SRC=../../src/integrator_whfast512.s
NSTEPS=${NSTEPS:-110000000000}
NSAMP=${NSAMP:-110}
NSEEDS=${NSEEDS:-80}

build_into(){
    cp "asm_variants/$1" "$SRC"
    make CC=gcc > /dev/null 2>&1
    test -x rebound || { echo "BUILD FAILED for $1"; exit 1; }
    mkdir -p "$2"; rm -f "$2"/*.txt
    cp rebound librebound.so "$2"/
}

build_into wf_fixed.s run_fixed
build_into wf_baseline.s run_baseline
cp asm_variants/wf_fixed.s "$SRC"

echo "launching $((2*(NSEEDS+1))) sims, nproc=$(nproc), NSTEPS=$NSTEPS"
for s in $(seq 0 $NSEEDS); do
    ( cd run_fixed    && LD_LIBRARY_PATH=. ./rebound $s $NSTEPS $NSAMP 1 17 1000 whfast512 6 > "seed$s.txt" 2>/dev/null ) &
    ( cd run_baseline && LD_LIBRARY_PATH=. ./rebound $s $NSTEPS $NSAMP 1 17 1000 whfast512 6 > "seed$s.txt" 2>/dev/null ) &
done
wait
echo "done: run_fixed=$(ls run_fixed/*.txt|wc -l) run_baseline=$(ls run_baseline/*.txt|wc -l)"
