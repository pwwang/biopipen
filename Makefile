test_targets := $(patsubst tests/test_%.py,%,$(wildcard tests/test_*.py))

test: $(test_targets)

%: tests/test_%.py
	@echo $<
	@# python $<

.PHONY: test
