# This Makefile compiles the shared dynamic library to access the IAS15 integrator.
export

librebound: 
	$(MAKE) -C src librebound
	@cp src/librebound.so .
	@if [ "$(MAKELEVEL)" -eq "0" ]; then echo "To compile the example problems, go to a subdirectory of examples/ and execute make there."; fi
	
all: librebound
