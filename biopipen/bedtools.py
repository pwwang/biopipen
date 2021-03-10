from .core.proc import Proc
from .core.args import args
from .core.defaults import SCRIPT_DIR

MODULE = 'bedtools'

class Coverage(Proc):
    """Run `bedtools coverage`

    See: https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html

    Input:
        afile: File fed to `-a` option
        bfile: File fed to `-b` option

    Output:
        outfile: The file reporting the coverage

    Args:
        bedtools: Path to bedtools
        params: Other params for `bedtools coverage`
    """
    lang = args.python
    input_keys = 'afile:file, bfile:file'
    output = 'outfile:file:{{in.bfile | stem}}-{{in.afile | stem}}.coverage'
    script = f'file://{SCRIPT_DIR}/{MODULE}/Coverage.py'
    args = {
        'bedtools': args.bedtools,
        'params': {}
    }
