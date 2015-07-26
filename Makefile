# This Makefile compiles the shared dynamic library to access the IAS15 integrator.
export OPT+= -std=c99 -Wpointer-arith -D_GNU_SOURCE -O3 -march=native

librebound: 
	$(MAKE) -C src 
	@ln -f -s src/librebound.so .
	@if [ "$(MAKELEVEL)" -eq "0" ]; then echo "To compile the example problems, go to a subdirectory of examples/ and execute make there."; fi
	
all: librebound

clean:
	$(MAKE) -C src clean
