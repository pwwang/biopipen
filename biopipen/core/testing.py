"""Provide utilities for testing."""
import tempfile
from functools import wraps
from pathlib import Path

from pipen import Pipen

TESTING_INDEX_INIT = 1
TESTING_PARENT_DIR = Path(__file__).parent.parent.parent.joinpath("tests", "running")
TESTING_PARENT_DIR.mkdir(parents=True, exist_ok=True)
TESTING_DIR = str(TESTING_PARENT_DIR.joinpath("biopipen-tests-%(index)s"))
RSCRIPT_DIR = TESTING_PARENT_DIR.joinpath("biopipen-tests-rscripts")
RSCRIPT_DIR.mkdir(exist_ok=True)


def _find_testing_index(new):
    """Find the next available testing index"""
    index = TESTING_INDEX_INIT
    while True:
        dir = TESTING_DIR % {"index": index}
        if not Path(dir).exists():
            if new:
                break
            else:
                return max(index - 1, TESTING_INDEX_INIT)
        index += 1
    return index


def _get_test_dirs(testfile, new):
    """Get the workdir and outdir for a test pipeline"""
    index = _find_testing_index(new)
    workdir = TESTING_DIR % {"index": index}
    procname = Path(testfile).parent.stem
    nsname = Path(testfile).parent.parent.stem
    name = f"{nsname}.{procname}"
    outdir = f"{workdir}/{nsname}/output"
    workdir = f"{workdir}/{procname}/pipen"
    Path(workdir).mkdir(parents=True, exist_ok=True)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    return name, workdir, outdir


def get_pipeline(testfile, loglevel="debug", enable_report=False, **kwargs):
    """Get a pipeline for a test file"""
    name, workdir, outdir = _get_test_dirs(testfile, False)
    report_plugin_prefix = "+" if enable_report else "-"
    plugins = kwargs.pop("plugins", [])
    if any("report" in p for p in plugins if isinstance(p, str)):
        raise ValueError(
            "Do not pass `report` plugin to `get_pipeline(plugins=[...])`, "
            "use `enable_report` instead."
        )
    plugins.append(f"{report_plugin_prefix}report")
    kws = {
        "name": name,
        "workdir": workdir,
        "outdir": outdir,
        "loglevel": loglevel,
        "plugins": plugins,
    }
    kws.update(kwargs)
    return Pipen(**kws)


def _run_rcode(rcode: str) -> str:
    """Run R code and return the output"""
    import hashlib
    import textwrap
    import subprocess as sp

    # Use sha256 of rcode to name the file
    rcode_hash = hashlib.sha256(rcode.encode()).hexdigest()
    script_file = RSCRIPT_DIR.joinpath(f"rcode-{rcode_hash}.R")
    script_file.write_text(rcode)
    p = sp.Popen(["Rscript", str(script_file)], stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        out = (
            f"R codefile:\n  {script_file}\n"
            f"Error:\n{textwrap.indent(err.decode(), '  ')}"
        )
        return out

    return out.decode().strip()


def r_test(mem: callable) -> callable:
    """A decorator to test R code"""
    @wraps(mem)
    def decorator(self, *args, **kwargs):
        rcode = mem(self, *args, **kwargs)
        source = getattr(self, "SOURCE_FILE", None)
        expect = (
            "expect <- function(expr, ...) {\n"
            "  if (!expr) {\n"
            "    msg <- lapply(\n"
            "      list(...),\n"
            "      function(x) { ifelse(is.null(x), 'NULL', x) }\n"
            "    )\n"
            "    stop(paste0(unlist(msg), collapse = ' '))\n"
            "  }\n"
            "}\n"
        )
        rcode = f"{expect}\n\n{rcode}\n\ncat('PASSED')\n"
        if source is not None:
            if not isinstance(source, (list, tuple)):
                source = [source]

            libs = "\n".join([f"suppressWarnings(source('{s}'))" for s in source])
            rcode = f'{libs}\n\n{rcode}'

        out = _run_rcode(rcode)
        self.assertEqual(
            out,
            "PASSED",
            "\n-----------------------------\n"
            f"{out}"
            "\n-----------------------------\n"
        )

    return decorator
