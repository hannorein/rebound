PROBLEMDIR=tests/$(shell basename `pwd`)
export OPENGL=0
export OPT=-O3 
export QUADRUPOLE

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

test:
	make all QUADRUPOLE=1
	mv nbody nbody-quad
	make all QUADRUPOLE=0
	mv nbody nbody-mono
	rm -f error.txt
	./nbody-quad 0.1
	./nbody-quad 0.2
	./nbody-quad 0.4
	./nbody-quad 0.6
	./nbody-quad 0.8
	./nbody-quad 1.0
	./nbody-quad 1.2
	./nbody-quad 1.4
	./nbody-quad 1.6
	./nbody-quad 1.8
	./nbody-quad 2.0
	mv error.txt error-quad.txt
	./nbody-mono 0.1
	./nbody-mono 0.2
	./nbody-mono 0.4
	./nbody-mono 0.6
	./nbody-mono 0.8
	./nbody-mono 1.0
	./nbody-mono 1.2
	./nbody-mono 1.4
	./nbody-mono 1.6
	./nbody-mono 1.8
	./nbody-mono 2.0
	cp error.txt error-mono.txt
	cat error-quad.txt >> error.txt
	gnuplot plot.plt

