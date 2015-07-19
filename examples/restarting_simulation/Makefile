PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
#export OPENGL=1
export OPT+= -g

all:
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/rebound .

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf rebound
