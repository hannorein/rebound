PROBLEMDIR=tests/$(shell basename `pwd`)
export OPENGL=1
export OPT=-O3

all:
	# Setup link to different modules
	ln -fs gravity_tree.c ../../src/gravity.c
	ln -fs boundaries_shear.c ../../src/boundaries.c
	ln -fs integrator_sei.c ../../src/integrator.c
	ln -fs collisions_sweep.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

test: all 
	rm -f error.txt
	./nbody 0.3 88.50 
	./nbody 0.5 19.10
	./nbody 0.6 11.05
	./nbody 0.65 8.69
	./nbody 0.675 8.0
	./nbody 0.7 6.97
	./nbody 0.75 3.14
	./nbody 0.8 6.01
	./nbody 0.82 6.01
	./nbody 0.85 6.04
	./nbody 0.9 6.0
	./nbody 0.95 6.04
	./nbody 1 6.04
