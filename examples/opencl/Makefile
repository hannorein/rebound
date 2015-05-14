PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export OPENGL=1
export OPT=-O3 -Wno-deprecated-declarations -g 
export LIB=-framework OpenCL
export QUADRUPOLE=0
export MPI=0

all:
	# Setup link to different modules
	ln -fs gravity_opencl.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/rebound .
direct:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/rebound .

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf rebound
