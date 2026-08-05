export OPENGL=0
export SERVER=0
include ../../src/Makefile.defs

.PHONY: all librebound clean

all: problem.c librebound
	$(CCPROBLEM)

librebound:
	$(MAKE) -C ../../src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)

clean:
	$(MAKE) -C ../../src/ clean
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)