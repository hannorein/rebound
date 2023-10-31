#!/bin/bash

source emsdk/emsdk_env.sh
READTHEDOCS_OUTPUT="${READTHEDOCS_OUTPUT:-.}"
OPENGL="${1:-1}"
OPTIMI="${2:-3}"
    
echo "Compiling C examples with emscripten."
echo "Output dir: $READTHEDOCS_OUTPUT"
echo "OPENGL: $OPENGL"
echo ""

for dir in examples/*/
do
    echo "Working on $dir ..."
    opengl_enabled=$(cat $dir/Makefile | grep -c "export OPENGL=1")
    mpi_enabled=$(cat $dir/Makefile | grep -c "export MPI=1")
    openmp_enabled=$(cat $dir/Makefile | grep -c "export OPENMP=1")
    if [ $opengl_enabled -eq $OPENGL ] && [ $mpi_enabled -eq 0 ] && [ $openmp_enabled -eq 0 ]; then
        mkdir -p $READTHEDOCS_OUTPUT/html/emscripten_c_$dir/
        echo "Compiling... "
        if [ $OPENGL -eq 1 ]; then
            emcc -O$OPTIMI -Isrc/ src/*.c $dir/problem.c -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY --shell-file docs/c_examples/shell_rebound.html -o $READTHEDOCS_OUTPUT/html/emscripten_c_$dir/index.html || exit 1
        else
            emcc -O$OPTIMI -Isrc/ src/*.c $dir/problem.c -sSTACK_SIZE=655360 -sASYNCIFY -o $READTHEDOCS_OUTPUT/html/emscripten_c_$dir/index.html || exit 1
        fi
        echo "Done. "
    else
        echo "Skipping."
    fi
    echo ""

done
