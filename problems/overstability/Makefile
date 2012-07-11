PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
#export CC=/Users/rein/localinstalls/src/gcc-lion/usr/local/bin/gcc
#export MPI=0
export OPENGL=1
export OPENMP=0
export OPT=-O3 -g 
export OPT+=-I/Users/rein/localinstalls/include
export LIB=
export LIB+=-L/Users/rein/localinstalls/lib -lfftw3_threads -lpthread 
export FFTW=1

all:
	# Setup link to different modules
	ln -fs gravity_none.c ../../src/gravity.c
	ln -fs boundaries_shear.c ../../src/boundaries.c
	ln -fs integrator_sei.c ../../src/integrator.c
	ln -fs collisions_sweep.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

doc: all
	$(MAKE) -C ../../src/ doc

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf nbody

