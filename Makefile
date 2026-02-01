SHELL=/bin/bash

NS_TARGETS := $(subst tests/test_,,$(wildcard tests/test_*))
PROC_TARGETS := $(subst tests/test_,,$(foreach ns,$(NS_TARGETS),$(wildcard tests/test_$(ns)/*)))

all: $(NS_TARGETS)

list:
	$(info Tests for namespaces: $(NS_TARGETS))
	$(info Tests for processes: $(PROC_TARGETS))

%: tests/test_%
	@echo "::group::Running tests for namespace: $@";                             \
	for procdir in $</*; do                                                       \
		echo "::group::├─ $@/$$(basename $$procdir)";                             \
		bash tests/conda/run_test.sh $$procdir VERBOSE=$(VERBOSE) FORCE=$(FORCE); \
		if [ $$? -ne 0 ]; then                                                    \
			should_fail=1;                                                        \
		fi;                                                                       \
		echo "::endgroup::";                                                      \
	done;                                                                         \
	echo "::endgroup::";                                                          \
	if [ -n "$$should_fail" ]; then                                               \
		exit 1;                                                                   \
	fi;

$(PROC_TARGETS): %: tests/test_%
	@bash tests/conda/run_test.sh $< VERBOSE=$(VERBOSE) FORCE=$(FORCE);

version:
	@if [ -z "$(word 2,$(MAKECMDGOALS))" ]; then \
		CURRENT_VERSION=$$(grep '^__version__' biopipen/__init__.py | sed 's/__version__ = "\(.*\)"/\1/'); \
		MAJOR=$$(echo $$CURRENT_VERSION | cut -d. -f1); \
		MINOR=$$(echo $$CURRENT_VERSION | cut -d. -f2); \
		PATCH=$$(echo $$CURRENT_VERSION | cut -d. -f3); \
		NEW_PATCH=$$((PATCH + 1)); \
		NEW_VERSION="$$MAJOR.$$MINOR.$$NEW_PATCH"; \
	else \
		NEW_VERSION="$(word 2,$(MAKECMDGOALS))"; \
	fi; \
	echo "Updating version to $$NEW_VERSION"; \
	sed -i "s/^version = .*/version = \"$$NEW_VERSION\"/" pyproject.toml; \
	sed -i "s/^__version__ = .*/__version__ = \"$$NEW_VERSION\"/" biopipen/__init__.py; \
	LAST_TAG=$$(git describe --tags --abbrev=0 2>/dev/null || echo ""); \
	if [ -z "$$LAST_TAG" ]; then \
		COMMITS=$$(git log --pretty=format:"- %s" HEAD); \
	else \
		COMMITS=$$(git log --pretty=format:"- %s" $$LAST_TAG..HEAD); \
	fi; \
	if [ -n "$$COMMITS" ]; then \
		printf "\n## %s\n\n%s\n\n" "$$NEW_VERSION" "$$COMMITS" | cat - <(tail -n +3 docs/CHANGELOG.md) > docs/CHANGELOG.md.tmp; \
		head -n 2 docs/CHANGELOG.md > docs/CHANGELOG.md.new; \
		cat docs/CHANGELOG.md.tmp >> docs/CHANGELOG.md.new; \
		mv docs/CHANGELOG.md.new docs/CHANGELOG.md; \
		rm -f docs/CHANGELOG.md.tmp; \
	fi; \
	echo "Version updated to $$NEW_VERSION";

# Catch-all rule to ignore version number argument
%:
	@:

.PHONY: all list version
