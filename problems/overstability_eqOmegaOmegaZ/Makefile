PROBLEMDIR=problems/$(shell basename `pwd`)
export CC=gcc
export OPENGL=0
export OPENMP=0
export OPT=-O3

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

run: all
	rm -rf out_tau*
	./nbody --tau=0.1
	./nbody --tau=0.2
	./nbody --tau=0.3
	./nbody --tau=0.4
	./nbody --tau=0.5
	./nbody --tau=0.6
	./nbody --tau=0.7
	./nbody --tau=0.8
	./nbody --tau=0.9
	./nbody --tau=1.0
	./nbody --tau=1.1
	./nbody --tau=1.2
	./nbody --tau=1.3
	./nbody --tau=1.4
	./nbody --tau=1.5
	./nbody --tau=1.6
	./nbody --tau=1.7
	./nbody --tau=1.8
	./nbody --tau=1.9
	./nbody --tau=2.0
	./nbody --tau=2.1
	./nbody --tau=2.2
	./nbody --tau=2.3
	./nbody --tau=2.4
	./nbody --tau=2.5
	./nbody --tau=2.6
	./nbody --tau=2.7
	./nbody --tau=2.8
	./nbody --tau=2.9
	./nbody --tau=3.0
	
