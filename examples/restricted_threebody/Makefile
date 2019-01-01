export OPENGL=1

# OpenMP
# On Mac OSX, we can use the CLANG compiler. But it requires some additional 
# flags (see Makefile.defs in src/ directory). You also need to install the 
# OpenMP library with homebrew:
#    brew install libomp
# Alternatively use a compiler which supports OpenMP out of the box (gcc) and
# uncomment the following line:
# export CC=gcc

# Ucomment the following lines to turn on OPENMP
#ifeq ($(shell $(CC) -v 2>&1 | grep -c "clang"), 1)
#export OPENMPCLANG=1
#else
#export OPENMP=1
#endif


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
