# This Makefile compiles the shared dynamic library librebound.so

librebound: 
	$(MAKE) -C src 
	@ln -f -s src/librebound.so .
	@if [ "$(MAKELEVEL)" -eq "0" ]; then echo "To compile the example problems, go to a subdirectory of examples/ and execute make there."; fi
	-cp librebound.so `python -c "import rebound; print(rebound.__libpath__)"`
	
all: librebound

clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean

.PHONY: doc
doc: 
	cd doc/doxygen && doxygen
	$(MAKE) -C doc html
		
