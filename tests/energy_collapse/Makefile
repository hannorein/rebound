PROBLEMDIR=tests/$(shell basename `pwd`)
export OPENGL=0
export OPT=-O3 

all:
	# Setup link to different modules
	ln -fs gravity_tree.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs integrator_leapfrog.c ../../src/integrator.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

test: all
	rm -f error.txt
	./nbody 0.2
	./nbody 0.4
	./nbody 0.6
	./nbody 0.8
