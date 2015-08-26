export OPENGL=0
include ../src/Makefile.defs

REBOUNDXFDIR=../../reboundxf/cext/

eclipse: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -g -I../src/ -I../../reboundxf/cext/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c $(REBOUNDXFDIR)libreboundxf.c $(REBOUNDXFDIR)xftools.c -L. -L ../../reboundxf -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

all: librebound libreboundxf
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -g -I../src/ -I../../reboundxf/cext/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -L ../../reboundxf -lreboundxf -lrebound $(LIB) -o rebound
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
