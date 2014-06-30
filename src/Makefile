OPT+= -std=c99 -Wpointer-arith
ifndef OS
	OS=$(shell uname)
endif
ifeq ($(OS), Linux)
	OPT+= -Wall  -g 
	LIB+= -lm
endif
ifeq ($(OS), Darwin)
	OPT+= -I/usr/local/include  -Wall -g
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
PREDEF+= -D$(shell basename `readlink integrator.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF+= -D$(shell basename `readlink boundaries.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF+= -D$(shell basename `readlink collisions.c` '.c' | tr '[a-z]' '[A-Z]')

PREDEFLIST= $(shell echo "$(PREDEF)" | sed 's/-D//g')

all:
	$(CC) $(OPT) $(PREDEF) main.c -c
	$(CC) $(OPT) $(PREDEF) problem.c -c
	$(CC) $(OPT) $(PREDEF) tree.c -c
	$(CC) $(OPT) $(PREDEF) particle.c -c
	$(CC) $(OPT) $(PREDEF) gravity.c -c
	$(CC) $(OPT) $(PREDEF) integrator.c -c
	$(CC) $(OPT) $(PREDEF) boundaries.c -c
	$(CC) $(OPT) $(PREDEF) input.c -c
	$(CC) $(OPT) $(PREDEF) output.c -c
	$(CC) $(OPT) $(PREDEF) collisions.c -c
	$(CC) $(OPT) $(PREDEF) collision_resolve.c -c
	$(CC) $(OPT) $(PREDEF) communication_mpi.c -c
	$(CC) $(OPT) $(PREDEF) zpr.c -c
	$(CC) $(OPT) $(PREDEF) display.c -c
	$(CC) $(OPT) $(PREDEF) tools.c -c
	$(CC) *.o -o nbody $(LIB)

clean:
	rm -vf *.o
	rm -vf nbody

doc:
	sed 's!PREDEFINED             =.*!PREDEFINED             = $(PREDEFLIST)!' ../doc/Doxyfile.input > ../doc/Doxyfile 
	cd ../doc/ && doxygen


