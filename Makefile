# This Makefile compiles the shared dynamic library to access the IAS15 integrator.

export CC=cc
export OPT=-march=native

librebound: 
	$(MAKE) -C src librebound
	@cp src/librebound.so .
	@echo "(To compile the example problems, go to a subdirectory of examples/ and execute make there.)"
	@echo ""        
	
all: librebound
