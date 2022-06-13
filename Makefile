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
		bash tests/conda/run_test.sh $$procdir $(VERBOSE); \
	done || exit 1;                                        \
	echo "::endgroup::"

$(PROC_TARGETS): %: tests/test_%
	@bash tests/conda/run_test.sh $< $(VERBOSE);

log:
	@lasttag=$$(git describe --tags --abbrev=0);                        \
	git --no-pager log --pretty=format:"%s" --reverse $$lasttag..HEAD | \
	while read line; do                                                 \
		echo $$line | tr ";" "\n" | while read log; do                  \
			log=$$(echo $$log | sed 's/^\s\+//');                       \
			if [[ "$$log" =~ ^[^:alnum:] ]]; then                       \
				echo "- $$log";                                         \
			fi;                                                         \
		done;                                                           \
	done;


.PHONY: all .list
