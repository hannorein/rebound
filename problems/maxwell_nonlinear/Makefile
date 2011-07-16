PROBLEMDIR=problems/maxwell_nonlinear
export OPENGL=1
export LIBPNG=1

all:
	# Setup link to different modules
	ln -fs ../$(PROBLEMDIR)/gravity_direct_cutoff.c ../../src/gravity.c
	ln -fs boundaries_sheary.c ../../src/boundaries.c
	ln -fs integrator_sei.c ../../src/integrator.c
	ln -fs collisions_direct.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .
