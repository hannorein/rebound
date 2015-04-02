# This Makefile compiles the shared dynamic library to access the IAS15 and MIKKOLA integrators.
CC=cc 
OPT=-std=c99 -O3 -DLIBREBOUND -DCOLLISIONS_NONE

all: 
	$(CC) $(OPT) -c -fPIC -Wall ../../src/integrator_ias15.c
	$(CC) $(OPT) -c -fPIC -Wall ../../src/integrator_mikkola.c
	$(CC) $(OPT) -c -fPIC -Wall ../../src/integrator_wh.c
	$(CC) $(OPT) -c -fPIC -Wall ../../src/integrator_leapfrog.c
	$(CC) $(OPT) -c -fPIC -Wall ../../src/gravity_direct.c
	$(CC) $(OPT) -c -fPIC -Wall ../../src/tools.c
	$(CC) $(OPT) -c -fPIC -Wall ../../src/particle.c
	$(CC) $(OPT) -c -fPIC -I../../src/ -Wall librebound.c
	$(CC) $(OPT) -shared -o librebound.so particle.o tools.o integrator_ias15.o integrator_mikkola.o integrator_wh.o integrator_leapfrog.o gravity_direct.o librebound.o 
	
	@echo ""	
	@echo "The shared library librebound.so has been created. Cleaning up"
	@echo ""	

	rm -f *.o

