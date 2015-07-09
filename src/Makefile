OPT+= -std=c99 -Wpointer-arith
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
	CC?=gcc
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

ifeq ($(PARTICLEIDS), 1)
	PREDEF+= -DPARTICLEIDS
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



# Create precompiler definitions from links
PREDEF+= -D$(shell basename `readlink gravity.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF+= -D$(shell basename `readlink boundaries.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF+= -D$(shell basename `readlink collisions.c` '.c' | tr '[a-z]' '[A-Z]')

PREDEFLIST= $(shell echo "$(PREDEF)" | sed 's/-D//g')

all:
	# Remove old object files
	rm -f *.o
	# Compile all source files 
	@if test -L integrator.c; then { echo "\n\nERROR! The file integrator.c seems to be a symbolic link. In this version of REBOUND you can choose the integrator in your problem.c file and no longer have to overwrite the integrator.c file. Please check the README file and update your Makefile accordingly. Exiting...\n\n"; exit 1;}; else exit 0; fi
	$(CC) $(OPT) $(PREDEF) main.c -c
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
	$(CC) $(OPT) $(PREDEF) boundaries.c -c
	$(CC) $(OPT) $(PREDEF) input.c -c
	$(CC) $(OPT) $(PREDEF) output.c -c
	$(CC) $(OPT) $(PREDEF) collisions.c -c
	$(CC) $(OPT) $(PREDEF) collision_resolve.c -c
	$(CC) $(OPT) $(PREDEF) communication_mpi.c -c
	$(CC) $(OPT) $(PREDEF) zpr.c -c
	$(CC) $(OPT) $(PREDEF) display.c -c
	$(CC) $(OPT) $(PREDEF) tools.c -c
	$(CC) *.o -o rebound $(LIB)

clean:
	rm -vf *.o
	rm -vf rebound


# This part of the Makefile compiles the shared dynamic library librebound.
librebound: OPT+=-DCOLLISIONS_NONE -DLIBREBOUND
librebound:
	$(CC) $(OPT) -c -fPIC integrator_ias15.c
	$(CC) $(OPT) -c -fPIC integrator_whfast.c
	$(CC) $(OPT) -c -fPIC integrator_leapfrog.c
	$(CC) $(OPT) -c -fPIC integrator_hybrid.c
	$(CC) $(OPT) -c -fPIC integrator_sei.c
	$(CC) $(OPT) -c -fPIC integrator_wh.c
	$(CC) $(OPT) -c -fPIC integrator.c
	$(CC) $(OPT) -c -fPIC gravity_direct.c
	$(CC) $(OPT) -c -fPIC tools.c
	$(CC) $(OPT) -c -fPIC particle.c
	$(CC) $(OPT) -c -fPIC librebound.c
	$(CC) $(OPT) -c -fPIC output.c
	$(CC) $(OPT) -c -fPIC input.c
	$(CC) $(OPT) -shared -o librebound.so particle.o tools.o integrator.o integrator_ias15.o integrator_whfast.o integrator_leapfrog.o integrator_sei.o integrator_hybrid.o integrator_wh.o gravity_direct.o librebound.o input.o output.o
	
	@echo ""        
	@echo "The shared library librebound.so has been created successfully."
	@echo ""        
	
	
