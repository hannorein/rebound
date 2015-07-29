# This Makefile compiles the shared dynamic library librebound.so

librebound: 
	$(MAKE) -C src 
	@ln -f -s src/librebound.so .
	@if [ "$(MAKELEVEL)" -eq "0" ]; then echo "To compile the example problems, go to a subdirectory of examples/ and execute make there."; fi
	
all: librebound

clean:
	$(MAKE) -C src clean
