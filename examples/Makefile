ifndef OS
	OS=$(shell uname)
endif
MPIEXAMPLES := shearing_sheet_mpi/. selfgravity_disc_mpi/. mpi_unittests/.
WHFAST512EXAMPLES := whfast512_unittests/. whfast512_solar_system/. whfast512_2planets/.
ifneq ($(OS), Windows_NT)
SUBDIRS := $(filter-out $(WHFAST512EXAMPLES) $(MPIEXAMPLES),$(wildcard */.))
else
SUBDIRS := $(filter-out openmp/. $(WHFAST512EXAMPLES) $(MPIEXAMPLES),$(wildcard */.))
endif
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

