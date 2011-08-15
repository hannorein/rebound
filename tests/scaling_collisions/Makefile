PROBLEMDIR=tests/$(shell basename `pwd`)
export OPENGL=0
export OPT=-O3

all:
	# Setup link to different modules
	ln -fs gravity_none.c ../../src/gravity.c
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

test: all
	./nbody 12	
	./nbody 25	
	./nbody 50	
	./nbody 100	
	./nbody 200	
	./nbody 400	
	./nbody 800	
	./nbody 1600	
	./nbody 3200	
	./nbody 5400	
	./nbody 12800	
