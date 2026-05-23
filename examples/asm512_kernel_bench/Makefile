export AVX512=1
export OPENGL=0
export SERVER=0
export OPT=-march=native -O3 -fmax-errors=5
include ../../src/Makefile.defs

.PHONY: all librebound
all: problem.c librebound
	@echo "Compiling $< ..."
	$(CCPROBLEM)
	@echo "Compilation successful. Run with: ./rebound [bench|chunks]"

librebound:
	@$(MAKE) -C ../../src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)

clean:
	@$(MAKE) -C ../../src/ clean
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)
