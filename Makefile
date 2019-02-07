SHELL=bash
PYTHON=python
PYTHON2=python2
PYTHON3=python3

.PHONY: testutils tu testprocs tp

tu: ./tests/utils/test*.py
	@wd=`pwd`;                                                                         \
	cd ./tests/utils;                                                                  \
	for tfile in test*.py; do                                                          \
		echo "";                                                                       \
		echo "Running test: $(PYTHON) $$tfile";                                        \
		echo '----------------------------------------------------------------------'; \
		$(PYTHON) $$tfile;                                                             \
		if [[ $$? -ne 0 ]]; then exit 1; fi;                                           \
	done;                                                                              \
	cd $$wd

tp: ./tests/procs/test*.py
	@wd=`pwd`;                                                                         \
	cd ./tests/procs;                                                                  \
	for tfile in test*.py; do                                                          \
		echo "";                                                                       \
		echo "Running test: $(PYTHON) $$tfile";                                        \
		echo '----------------------------------------------------------------------'; \
		$(PYTHON) $$tfile;                                                             \
		if [[ $$? -ne 0 ]]; then exit 1; fi;                                           \
	done;                                                                              \
	cd $$wd

testutils: tu
testprocs: tp