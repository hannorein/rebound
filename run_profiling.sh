#!/bin/bash
cd examples/whfast512_solar_system_jac
export OPT="-march=native -O3 -g -fno-omit-frame-pointer"
make clean
make

cd ../..

echo "--- IPC Profiling ---" > profiling_output.txt
perf stat -e instructions,cycles python benchmark_driver.py --filter "WHFast512 Jacobi (Baseline)" > /dev/null 2>> profiling_output.txt

echo "--- Cache Profiling ---" >> profiling_output.txt
perf stat -e L1-dcache-loads,L1-dcache-load-misses python benchmark_driver.py --filter "WHFast512 Jacobi (Baseline)" > /dev/null 2>> profiling_output.txt

echo "--- Perf Record ---" >> profiling_output.txt
perf record --call-graph fp python benchmark_driver.py --filter "WHFast512 Jacobi (Baseline)" > /dev/null 2>> profiling_output.txt

perf annotate --stdio > perf_annotation.txt
