MPIEXAMPLES := shearing_sheet_mpi/. selfgravity_disc_mpi/. mpi_unittests/.
SUBDIRS := $(filter-out $(MPIEXAMPLES),$(wildcard */.))
MPISUBDIRS := $(filter $(MPIEXAMPLES),$(wildcard */.))

all: $(SUBDIRS)
$(SUBDIRS):
	@echo "Trying to compile example $(subst /.,,$@):"
	$(MAKE) -s -C $@ clean 
	$(MAKE) OPENGL=0 -j -s -C $@ 
	@echo "\033[1A\033[55CSuccess."

mpi: $(MPISUBDIRS)
$(MPISUBDIRS):
	@echo "Trying to compile MPI example $(subst /.,,$@):"
	$(MAKE) -s -C $@ clean 
	$(MAKE) OPENGL=0 -j -s -C $@ 
	@echo "\033[1A\033[55CSuccess."

.PHONY: all $(SUBDIRS)

.PHONY: mpi $(MPISUBDIRS)

