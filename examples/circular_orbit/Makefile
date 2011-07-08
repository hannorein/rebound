PROBLEMDIR=examples/circular_orbit

all:
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_periodic.c ../../src/boundaries.c
	ln -fs integrator_leapfrog.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	$(MAKE) -C ../../src/
	cp ../../src/nbody .
