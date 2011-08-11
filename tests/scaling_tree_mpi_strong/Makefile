PROBLEMDIR=tests/$(shell basename `pwd`)
export OPT=-O3
export QUADRUPOLE=0
export OPENGL=0
export MPI=1
export OPENMP=0

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

doc: all
	$(MAKE) -C ../../src/ doc

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf nbody

test:
	make all 
	rm -vf scaling.txt
	./nbody 10.0
	./nbody 20.0
	./nbody 40.0
	./nbody 80.0
	./nbody 160.0
	./nbody 320.0
	./nbody 640.0
	gnuplot plot.plt
