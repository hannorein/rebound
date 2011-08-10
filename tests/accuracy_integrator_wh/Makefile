PROBLEMDIR=tests/$(shell basename `pwd`)
export OPENGL=0

all:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs integrator_wh.c ../../src/integrator.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

test: all
	rm -f error.txt
	./nbody 0.0
	./nbody 0.1
	./nbody 0.2
	./nbody 0.3
	./nbody 0.4
	./nbody 0.5
	./nbody 0.6
	./nbody 0.7
	./nbody 0.8
	./nbody 0.9
	./nbody 0.99

