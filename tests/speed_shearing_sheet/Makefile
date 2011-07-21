PROBLEMDIR=tests/$(shell basename `pwd`)
export OPENGL=0
export OPT=-O3 
export QUADRUPOLE=1

all:
	# Setup link to different modules
	ln -fs gravity_tree.c ../../src/gravity.c
	ln -fs boundaries_shear.c ../../src/boundaries.c
	ln -fs integrator_sei.c ../../src/integrator.c
	ln -fs collisions_tree.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

test: all
	rm -f error.txt
	./nbody 0.5
	./nbody 1
	./nbody 2
	./nbody 4
