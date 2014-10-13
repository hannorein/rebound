PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export OPENGL=1
export OPT=-O3
export CC=gcc
export FORTRAN=gfortran
export LIB=../$(PROBLEMDIR)/additional_fortran_forces.o -L/usr/local/gfortran/lib/ -lgfortran

all:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs integrator_ias15.c ../../src/integrator.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c

	# Compile
	$(FORTRAN) -ffree-form -c additional_fortran_forces.f90
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .
