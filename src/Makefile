include Makefile.defs

SOURCES=rebound.c tree.c particle.c gravity.c integrator.c integrator_whfast.c integrator_whfast512.c integrator_saba.c integrator_ias15.c integrator_sei.c integrator_bs.c integrator_leapfrog.c integrator_mercurius.c integrator_trace.c integrator_eos.c boundary.c input.c binarydiff.c output.c collision.c communication_mpi.c display.c tools.c rotations.c derivatives.c simulationarchive.c glad.c integrator_janus.c transformations.c fmemopen.c server.c
OBJECTS=$(SOURCES:.c=.$(OBJFILEEXT))
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) $(LIBREBOUND)

ifneq ($(OS),Windows_NT) 
# Linux, MacOS
%.$(OBJFILEEXT): %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -c $(OPT) $(PREDEF) -o $@ $<

$(LIBREBOUND): $(OBJECTS)
	@echo ""        
	@echo "Linking shared library $@ ..."
	$(CC) $(OPT) -shared $(OBJECTS) $(LIB) -o $@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
else
# Windows
%.$(OBJFILEEXT): %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -c $(OPT) $(PREDEF) /Fo: $@ $<

$(LIBREBOUND): $(OBJECTS)
	@echo ""        
	@echo "Linking shared library $@ ..."
	$(CC) /D_USRDLL /D_WINDLL $(OBJECTS) /link /DLL /OUT:$@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
endif
	
	
clean:
	@echo "Removing object files *.$(OBJFILEEXT) ..."
	@-$(RM) *.$(OBJFILEEXT)
	@echo "Removing shared library $(LIBREBOUND) ..."
	@-$(RM) *.so
	@echo "Removing coverage files ..."
	@-$(RM) *.gcda
	@-$(RM) *.gcno
	
