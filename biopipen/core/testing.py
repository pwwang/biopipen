"""Provide utilities for testing."""

import tempfile
from pathlib import Path

from pipen import Pipen

TESTING_INDEX_INIT = 1
TESTING_PARENT_DIR = tempfile.gettempdir()
TESTING_DIR = f"{TESTING_PARENT_DIR}/biopipen-tests-%(index)s"


def _find_testing_index():
    """Find the next available testing index"""
    index = TESTING_INDEX_INIT
    while True:
        dir = TESTING_DIR % {"index": index}
        if not Path(dir).exists():
            break
        index += 1
    return index


def _get_test_dirs(testfile):
    """Get the workdir and outdir for a test pipeline"""
    index = _find_testing_index()
    workdir = TESTING_DIR % {"index": index}
    procname = Path(testfile).parent.stem
    nsname = Path(testfile).parent.parent.stem
    name = f"{nsname}.{procname}"
    outdir = f"{workdir}/{nsname}/output"
    workdir = f"{workdir}/{procname}/pipen"
    Path(workdir).mkdir(parents=True, exist_ok=True)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    return name, workdir, outdir


def get_pipeline(testfile, loglevel="debug", **kwargs):
    """Get a pipeline for a test file"""
    name, workdir, outdir = _get_test_dirs(testfile)
    kws = {
        "name": name,
        "workdir": workdir,
        "outdir": outdir,
        "loglevel": loglevel,
    }
    kws.update(kwargs)
    pipen = Pipen(**kws)
    return pipen
