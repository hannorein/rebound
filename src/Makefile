OPT+= -std=c99 -Wpointer-arith -D_GNU_SOURCE
ifndef OS
	OS=$(shell uname)
endif
ifeq ($(OS), Linux)
	OPT+= -Wall -g
	LIB+= -lm
endif
ifeq ($(OS), Darwin)
	OPT+= -I/usr/local/include -Wall -Wno-deprecated -g
	PREDEF+= -D_APPLE
	LIB+= -L/usr/local/lib
endif

ifeq ($(MPI), 1)
	CC?=mpicc
	PREDEF+= -DMPI
	OPENGL=0
else
	CC?=cc
endif

ifeq ($(FFTW), 1)
	PREDEF+= -DFFTW
	LIB+= -lfftw3
endif

ifeq ($(OPENGL), 1)
PREDEF+= -DOPENGL
ifeq ($(OS), Darwin)
	LIB+= -framework OpenGL -framework GLUT
else
	LIB+= -lGL -lGLU -lglut 
endif

ifeq ($(LIBPNG), 1)
	PREDEF+= -DLIBPNG
	LIB+= -lpng
endif
endif

ifeq ($(QUADRUPOLE), 1)
	PREDEF+= -DQUADRUPOLE
endif

ifeq ($(PROFILING), 1)
	PREDEF+= -DPROFILING
endif

ifeq ($(OPENMP), 1)
	PREDEF+= -DOPENMP
ifeq ($(CC), icc)
	OPT+= -openmp
	LIB+= -openmp
else
	OPT+= -fopenmp
	LIB+= -fopenmp
endif
else
	OPT+= -Wno-unknown-pragmas
endif



PREDEFLIST= $(shell echo "$(PREDEF)" | sed 's/-D//g')

all:
	# Remove old object files
	rm -f *.o
	# Compile all source files 
	@if test -L integrator.c; then { echo "\n\nERROR! The file integrator.c seems to be a symbolic link. In this version of REBOUND you can choose the integrator in your problem.c file and no longer have to overwrite the integrator.c file. Please check the README file and update your Makefile accordingly. Exiting...\n\n"; exit 1;}; else exit 0; fi
	$(CC) $(OPT) $(PREDEF) rebound.c -c
	$(CC) $(OPT) $(PREDEF) problem.c -c
	$(CC) $(OPT) $(PREDEF) tree.c -c
	$(CC) $(OPT) $(PREDEF) particle.c -c
	$(CC) $(OPT) $(PREDEF) gravity.c -c
	$(CC) $(OPT) $(PREDEF) integrator.c -c
	$(CC) $(OPT) $(PREDEF) integrator_whfast.c -c
	$(CC) $(OPT) $(PREDEF) integrator_ias15.c -c
	$(CC) $(OPT) $(PREDEF) integrator_sei.c -c
	$(CC) $(OPT) $(PREDEF) integrator_wh.c -c
	$(CC) $(OPT) $(PREDEF) integrator_leapfrog.c -c
	$(CC) $(OPT) $(PREDEF) integrator_hybrid.c -c
	$(CC) $(OPT) $(PREDEF) boundary.c -c
	$(CC) $(OPT) $(PREDEF) input.c -c
	$(CC) $(OPT) $(PREDEF) output.c -c
	$(CC) $(OPT) $(PREDEF) collision.c -c
	$(CC) $(OPT) $(PREDEF) communication_mpi.c -c
	$(CC) $(OPT) $(PREDEF) zpr.c -c
	$(CC) $(OPT) $(PREDEF) display.c -c
	$(CC) $(OPT) $(PREDEF) tools.c -c
	$(CC) *.o -o rebound $(LIB)

clean:
	rm -vf *.o
	rm -vf rebound


# This part of the Makefile compiles the shared dynamic library librebound.
librebound: OPT+= -fPIC -DLIBREBOUND

librebound:
	# Remove old object files
	rm -f *.o
	$(CC) $(OPT) $(PREDEF) rebound.c -c
	$(CC) $(OPT) $(PREDEF) tree.c -c
	$(CC) $(OPT) $(PREDEF) particle.c -c
	$(CC) $(OPT) $(PREDEF) gravity.c -c
	$(CC) $(OPT) $(PREDEF) integrator.c -c
	$(CC) $(OPT) $(PREDEF) integrator_whfast.c -c
	$(CC) $(OPT) $(PREDEF) integrator_ias15.c -c
	$(CC) $(OPT) $(PREDEF) integrator_sei.c -c
	$(CC) $(OPT) $(PREDEF) integrator_wh.c -c
	$(CC) $(OPT) $(PREDEF) integrator_leapfrog.c -c
	$(CC) $(OPT) $(PREDEF) integrator_hybrid.c -c
	$(CC) $(OPT) $(PREDEF) boundary.c -c
	$(CC) $(OPT) $(PREDEF) input.c -c
	$(CC) $(OPT) $(PREDEF) output.c -c
	$(CC) $(OPT) $(PREDEF) collision.c -c
	$(CC) $(OPT) $(PREDEF) communication_mpi.c -c
	$(CC) $(OPT) $(PREDEF) zpr.c -c
	$(CC) $(OPT) $(PREDEF) display.c -c
	$(CC) $(OPT) $(PREDEF) tools.c -c
	$(CC) $(OPT) -shared *.o -o librebound.so $(LIB)
	
	@echo ""        
	@echo "The shared library librebound.so has been created successfully."
	@echo ""        
	
	
