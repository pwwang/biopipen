SHELL=/bin/bash

NS_TARGETS := $(subst tests/test_,,$(wildcard tests/test_*))
PROC_TARGETS := $(subst tests/test_,,$(foreach ns,$(NS_TARGETS),$(wildcard tests/test_$(ns)/*)))

all: $(NS_TARGETS)

.list:
	$(info Tests for namespaces: $(NS_TARGETS))
	$(info Tests for processes: $(PROC_TARGETS))

%: tests/test_%
	@echo "::group::Running tests for namespace: $@";      \
	for procdir in $</*; do                                \
		bash tests/conda/run_test.sh $$procdir VERBOSE=$(VERBOSE) FORCE=$(FORCE); \
	done || exit 1;                                        \
	echo "::endgroup::"

$(PROC_TARGETS): %: tests/test_%
	@bash tests/conda/run_test.sh $< VERBOSE=$(VERBOSE) FORCE=$(FORCE);


.PHONY: all .list
