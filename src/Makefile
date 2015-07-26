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

OPT+= -fPIC -DLIBREBOUND

SOURCES=rebound.c tree.c particle.c gravity.c integrator.c integrator_whfast.c integrator_ias15.c integrator_sei.c integrator_wh.c integrator_leapfrog.c integrator_hybrid.c boundary.c input.c output.c collision.c communication_mpi.c zpr.c display.c tools.c 
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) librebound.so 

%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -c $(OPT) $(PREDEF) -o $@ $<

librebound.so: $(OBJECTS)
	@echo ""        
	@echo "Linking share library $@ ..."
	$(CC) $(OPT) -shared $(OBJECTS) $(LIB) -o $@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
	
	
clean:
	@echo "Removing object files *.o ..."
	@-rm -f *.o
	@echo "Removing shared library librebound.so ..."
	@-rm -f *.so
	
