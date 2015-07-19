SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)
$(SUBDIRS):
	@echo "Trying to compile example in $@:"
	@$(MAKE) -s -C $@ 
	@echo "Success."
	@echo ""
	@echo ""

.PHONY: all $(SUBDIRS)

