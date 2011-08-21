OPT+= -std=c99
ifeq ($(OS), Linux)
	OPT+= -I/usr/include -Wall  -g 
	LIB+= -L/usr/lib -lm
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

ifeq ($(OPENGL), 1)
ifeq ($(OS), Linux)
	PREDEF+= -DOPENGL
	LIB+= -lglut 
endif
ifeq ($(OS), Darwin)
	PREDEF+= -DOPENGL
	LIB+= -framework OPENGL -framework GLUT
endif

ifeq ($(LIBPNG), 1)
	PREDEF+= -DLIBPNG
	LIB+= -lpng
endif
endif

ifeq ($(QUADRUPOLE), 1)
	PREDEF+= -DQUADRUPOLE
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
	$(CC) $(LIB) *.o -o nbody

clean:
	rm -vf *.o
	rm -vf nbody

doc:
	sed 's!PREDEFINED             =.*!PREDEFINED             = $(PREDEFLIST)!' ../doc/Doxyfile.input > ../doc/Doxyfile 
	cd ../doc/ && doxygen


