include Makefile.defs

OPT+= -fPIC -DLIBREBOUND

SOURCES=rebound.c tree.c particle.c gravity.c integrator.c integrator_whfast.c integrator_whfasthelio.c integrator_ias15.c integrator_sei.c integrator_leapfrog.c integrator_hermes.c boundary.c input.c output.c collision.c communication_mpi.c display.c tools.c derivatives.c simulationarchive.c glad.c integrator_janus.c transformations.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) librebound.so 

%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -c $(OPT) $(PREDEF) -o $@ $<

librebound.so: $(OBJECTS)
	@echo ""        
	@echo "Linking share library $@ ..."
	$(CC) $(OPT) -shared $(OBJECTS) $(LIB) -o $@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
	
	
clean:
	@echo "Removing object files *.o ..."
	@-rm -f *.o
	@echo "Removing shared library librebound.so ..."
	@-rm -f *.so
	@echo "Removing coverage files ..."
	@-rm -f *.gcda
	@-rm -f *.gcno
	
