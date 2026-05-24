export OPENGL=0# Set this to 1 to enable OpenGL
export SERVER=1# Set this to 1 to enable the visualization web server
include ../../../../src/Makefile.defs

ifneq ($(OS), Windows_NT)
LINKORCOPYLIBREBOUND=ln -s -f ../../../../src/$(LIBREBOUND) .
CCPROBLEM=$(CC) -I../../../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $< -L. -lrebound $(LIB) -o $(EXEREBOUND)
else
LINKORCOPYLIBREBOUND=copy ..\..\..\..\src\librebound.* .
CCPROBLEM=$(CC) /I"..\..\..\..\src" $(OPT) $(PREDEF) $< librebound.lib /Fe: $(EXEREBOUND)
endif

.PHONY: all librebound
all: problem.c librebound
	@echo "Compiling $< ..."
	$(CCPROBLEM)
	@echo ""
	@echo "Compilation successful. To run REBOUND, execute the file '$(EXEREBOUND)'."
	@echo ""

librebound:
	@echo "Compiling shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../../../src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)
	@echo ""

clean:
	@echo "Cleaning up shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../../../src/ clean
	@echo "Cleaning up local directory ..."
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)
