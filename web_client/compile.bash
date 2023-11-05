#!/bin/bash
# To be called from main rebound directory
source emsdk/emsdk_env.sh
          
emcc -O3 -Isrc/ src/*.c web_client/problem.c -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sFETCH -sSINGLE_FILE --shell-file c_examples/shell_rebound.html -o web_client/rebound.html
