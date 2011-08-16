PROBLEMDIR=tests/$(shell basename `pwd`)
export OPENGL=0
export OPT=-O3

all:
	# Setup link to different modules
	ln -fs gravity_none.c ../../src/gravity.c
	ln -fs boundaries_shear.c ../../src/boundaries.c
	ln -fs integrator_sei.c ../../src/integrator.c
	ln -fs collisions_direct.c ../../src/collisions.c
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
	./nbody 1	
	./nbody 2	
	./nbody 4	
	./nbody 8	
	./nbody 16	
	./nbody 32	
	./nbody 64	
	./nbody 128	
	./nbody 256	
	./nbody 512	
	./nbody 1024	
