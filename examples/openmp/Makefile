# Turninng on OpenMP
# On Mac OSX, we can use the CLANG compiler. But it requires some additional
# flags (see Makefile.defs in src/ directory). You also need to install the
# OpenMP library with homebrew:
#    brew install libomp
# Alternatively use a compiler which supports OpenMP out of the box (gcc) and
# uncomment the following line:
#export CC=gcc

ifeq ($(shell $(CC) -v 2>&1 | grep -c "clang"), 1)
export OPENMPCLANG=1
else
export OPENMP=1
endif


export OPENGL=0  	# Set this to 1 to enable OpenGL
export SERVER=1		# Set this to 1 to enable the visualization web server
include ../../src/Makefile.defs

# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
.PHONY: all librebound
all: problem.c librebound
	@echo "Compiling $< ..."
	$(CCPROBLEM)
	@echo ""
	@echo "Compilation successful. To run REBOUND, execute the file '$(EXEREBOUND)'."
	@echo ""

librebound:
	@echo "Compiling shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)
	@echo ""

clean:
	@echo "Cleaning up shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)
