SHELL=/bin/bash

NS_TARGETS := $(patsubst tests/test_%,%,$(wildcard tests/test_*))


all: $(NS_TARGETS)

.list:
	@echo "Tests for namespaces: $(NS_TARGETS)"

%: tests/test_%
	@echo "Running tests for namespace: $@";         \
	for procdir in $</*; do                          \
		bash tests/conda/run_test.sh $$procdir $(VERBOSE); \
	done

$(NS_TARGETS).%: tests/test_$(NS_TARGETS)/%
	@bash tests/conda/run_test.sh $< $(VERBOSE)

.PHONY: all .list
