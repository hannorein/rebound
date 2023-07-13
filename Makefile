# This Makefile compiles the shared dynamic library librebound.so

librebound: 
	$(MAKE) -C src 
	@ln -f -s src/librebound.so .
	@if [ "$(MAKELEVEL)" -eq "0" ]; then echo "To compile the example problems, go to a subdirectory of examples/ and execute make there."; fi

.PHONY: pythoncopy
pythoncopy:
	-cp librebound.so `python -c "import rebound; print(rebound.__libpath__)"`
	
all: librebound pythoncopy

clean:
	$(MAKE) -C src clean
