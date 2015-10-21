export OPENGL=0
include ../src/Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -g -I../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../src/
	@-rm -f librebound.so
	@ln -s ../src/librebound.so .

libreboundxf: 
	@echo "Compiling shared library libreboundxf.so ..."
	$(MAKE) -C ../../reboundxf/cext/
	@-rm -f libreboundxf.so
	@ln -s ../../reboundxf/cext/libreboundxf.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	@-rm -f libreboundxf.so
	$(MAKE) -C ../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
