# turn on OpenMP
export OPENMP=1


# On Mac OSX, we'll have to use the GNU compilers.
# You can install them with homebrew. 
OS=$(shell uname)
	ifeq ($(OS), Darwin)
	# We'll use specifically the GNU compiler version 7 
	export CC=gcc-7
	# Use the OSX linker
	export OPT=-Wa,-q
endif

# Include the other definitions from the default makefile
include ../../src/Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../../src/
	@-rm -f librebound.so
	@ln -s ../../src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
