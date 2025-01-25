#!/bin/bash
if [ -f documentation_lock.txt ]; then
    echo "Documentation update already in progress. Exiting."
    exit
fi
touch documentation_lock.txt

git pull

repository_head=` git rev-parse HEAD`
documentation_head=$(<documentation_head.txt)

echo "Repository head:    $repository_head"
echo "Documentation head: $documentation_head"


if [[ "$repository_head" == "$documentation_head" ]]; then
    echo "Documentation is up-to-date."
    rm documentation_lock.txt
    exit
fi

echo "Updating documentation."

## UPDATE START

mkdocs build  
../emsdk/emsdk activate latest
source ../emsdk/emsdk_env.sh

OPTIMI="${1:-3}"
    
echo "Compiling C examples with emscripten."
echo ""
mkdir -p site/emscripten_c_examples/

for dir in examples/*/
do
    echo "Working on $dir ..."
    mpi_enabled=$(cat $dir/Makefile | grep -c "export MPI=1")
    openmp_enabled=$(cat $dir/Makefile | grep -c "export OPENMP=1")
    server_used=$(cat $dir/problem.c | grep -c "reb_simulation_start_server")
    if [ $mpi_enabled -eq 0 ] && [ $openmp_enabled -eq 0 ]; then
        mkdir -p site/emscripten_c_$dir/
        echo "Compiling... "
        if [ $server_used -eq 0 ]; then
            emcc -O$OPTIMI -Isrc/ src/*.c $dir/problem.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -s -sASYNCIFY -sALLOW_MEMORY_GROWTH -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file web_client/shell_rebound_console.html -o site/emscripten_c_$dir/index.html || exit 1
        else
            emcc -O$OPTIMI -Isrc/ src/*.c $dir/problem.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file web_client/shell_rebound_webgl.html -o site/emscripten_c_$dir/index.html || exit 1
        fi
        echo "Output written to site/emscripten_c_$dir/index.html" 
    else
        echo "Skipping."
    fi
    echo ""

done

## UPDATE DONE
echo "Documentaiton updated. Releasing lock."
echo "$repository_head" > documentation_head.txt
rm documentation_lock.txt

exit
