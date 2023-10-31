#!/bin/bash

source emsdk/emsdk_env.sh

for dir in examples/*
do
    echo $dir
    opengl_enabled=$(cat $dir/Makefile | grep -c "export OPENGL=1")
    mpi_enabled=$(cat $dir/Makefile | grep -c "export MPI=1")
    if [ $opengl_enabled -eq 1 ] && [ $mpi_enabled -eq 0 ]; then
        echo "Creating dir... "
        echo $READTHEDOCS_OUTPUT/html/emscripten_c_$dir/
        mkdir -p $READTHEDOCS_OUTPUT/html/emscripten_c_$dir/
        emcc -O3 -Isrc/ src/*.c $dir/problem.c -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sPTHREAD_POOL_SIZE=2 -s WASM_BIGINT -lwebsocket.js -sASYNCIFY --shell-file docs/c_examples/shell_rebound.html -o $READTHEDOCS_OUTPUT/html/emscripten_c_$dir/index.html
        echo "Compiling... "
    else
        echo "Skipping."
    fi
    echo ""

done
