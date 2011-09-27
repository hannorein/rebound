PROBLEMDIR=examples/$(shell basename `pwd`)
export OPT=-O3
export QUADRUPOLE=0
export OPENGL=1
export MPI=0
export FFTW=1
export OPENMP=0

all:
	# Setup link to different modules
	ln -fs gravity_fft.c ../../src/gravity.c
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
