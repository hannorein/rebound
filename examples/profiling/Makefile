export OPENGL=1
export PROFILING=1
OPT+= -std=c99 -Wpointer-arith -D_GNU_SOURCE -O3 -march=native -lm

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	cc -I../../src/ -L. -Wl,-rpath,./ $(OPT) problem.c -lrebound -o rebound
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
