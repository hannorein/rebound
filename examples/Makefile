SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)
$(SUBDIRS):
	@echo "Trying to compile example $(subst /.,,$@):"
	$(MAKE) -s -C $@ clean 
	$(MAKE) -j -s -C $@ 
	@echo "\033[1A\033[55CSuccess."

.PHONY: all $(SUBDIRS)

