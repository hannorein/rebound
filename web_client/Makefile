all:
	emcc -O3 -I../src/ ../src/*.c problem.c -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sGL_ENABLE_GET_PROC_ADDRESS=1 -sEXPORTED_RUNTIME_METHODS="callMain" -sFETCH -sSINGLE_FILE --shell-file shell_rebound.html -o rebound.html

clean:
	rm rebound.html
