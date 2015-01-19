# This Makefile compiles the shared dynamic library to access the IAS15 integrator.
CC=cc 
OPT=-std=c99 -O3 -DCOLLISIONS_NONE -DINTEGRATOR_IAS15

all: 
	$(CC) $(OPT) -c -fPIC -Wall ../src/integrator_ias15.c
	$(CC) $(OPT) -c -fPIC -Wall ../src/gravity_direct.c
	$(CC) $(OPT) -c -fPIC -Wall ../src/tools.c
	$(CC) $(OPT) -c -fPIC -Wall ../src/particle.c
	$(CC) $(OPT) -c -fPIC -I../src/ -Wall libias15.c
	$(CC) $(OPT) -shared -o libias15.so particle.o tools.o integrator_ias15.o gravity_direct.o libias15.o

	@echo ""	
	@echo "The shared library libias15.so has been created. Cleaning up"
	@echo ""	

	rm -f *.o

