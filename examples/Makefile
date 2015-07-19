SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)
$(SUBDIRS):
	@echo "Entering $@"
	$(MAKE) -C $@

.PHONY: all $(SUBDIRS)

