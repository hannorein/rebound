export OPENGL=0
export SERVER=0
export MPI=1
export CC=mpicc
include ../../src/Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound:
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../../src/
	@-rm -f librebound.so
	@ln -s ../../src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound

rebound_webgl.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL enabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_webgl.html -o rebound_webgl.html

rebound_console.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL disabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_console.html -o rebound_console.html
