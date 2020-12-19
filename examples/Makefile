SUBDIRS := $(filter-out shearing_sheet_mpi/. selfgravity_disc_mpi/.,$(wildcard */.))

all: $(SUBDIRS)
$(SUBDIRS):
	@echo "Trying to compile example $(subst /.,,$@):"
	$(MAKE) -s -C $@ clean 
	$(MAKE) OPENGL=0 -j -s -C $@ 
	@echo "\033[1A\033[55CSuccess."

.PHONY: all $(SUBDIRS)

