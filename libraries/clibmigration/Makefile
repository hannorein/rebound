# This Makefile compiles the shared dynamic library to access the IAS15 and MIKKOLA integrators.
CC=cc 
OPT=-std=c99 -O3 -DLIBREBOUND -DCOLLISIONS_NONE -DINTEGRATOR_NOTWH

all: 
	$(CC) $(OPT) -c -fPIC -I../clibrebound/ -I../../src/ -Wall libmigration.c
	$(CC) $(OPT) -L../clibrebound -shared -o libmigration.so libmigration.o -lrebound
	
	@echo ""	
	@echo "The shared library libmigration.so has been created. Cleaning up"
	@echo ""	

	rm -f *.o

