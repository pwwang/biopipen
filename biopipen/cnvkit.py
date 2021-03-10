from diot import Diot

from .core.proc import Proc
from .core.args import args
from .core.defaults import SCRIPT_DIR

MODULE = 'cnvkit'

class CNVkitBatch(Proc):
    """Run `cnvkit.py batch`

    See: https://cnvkit.readthedocs.io/en/stable/pipeline.html

    Input:
        tumors: Tumor bam files
        normals: Normal bam files. If not provided, a flat reference will
            be used

    Output:
        outdir: The output directory

    Args:
        cnvkit: Path to cnvkit.py
        access: The access file. If not provided, will use `cnvkit.py access`
            to generate
            If an integer is provided, will be used as `MIN_GAP_SIZE` for
            `cnvkit.py access`
        annotate: gene annotation file. (e.g http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz)
        params: Other params for `cnvkit.py batch`
    """
    lang = args.python
    input_keys = 'tumors:files, normals:files'
    output = 'outdir:file:{{in.tumors[0] | stem | append: ".cnvkit"}}'
    script = f'file://{SCRIPT_DIR}/{MODULE}/CNVkitBatch.py'
    args = Diot(
        cnvkit=args.cnvkit,
        ncores=1,
        ref=args.ref,
        params={
            'access': None,
            'targets': args.exbaits,
            'annotate': args.refflat,
            'diagram': True,
            'scatter': True,
        }
    )
