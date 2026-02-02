#!/bin/bash
set -e
cd examples/whfast512_solar_system_jac
make clean > /dev/null

export OPT="-march=native -O3 -DPROF"
# export OPT="-march=native -O3"

make
./rebound