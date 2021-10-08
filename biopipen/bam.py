from diot import Diot

from .core.proc import Proc
from .core.args import args
from .core.defaults import SCRIPT_DIR

MODULE = 'bam'

class PatternCNV(Proc):
    """Run `PatternCNV`

    Input:
        tumors: Tumor bam files
        normals: Normal bam files. If not provided, will treated as all germline
            calling

    Output:
        outdir: The output directory

    Args:
        patterncnv: Path to `patternCNV_wrapper.sh`
        samtools: Path to `samtools`
        bedtools: Path to `bedtools`
        r: Path to `R`
        params: Other params for `patternCNV_wrapper.sh`
    """
    lang = args.python
    input_keys = 'tumors:files, normals:files'
    output = 'outdir:file:{{in.tumors[0] | stem | append: ".patterncnv"}}'
    script = f'file://{SCRIPT_DIR}/{MODULE}/PatternCNV.py'
    args = Diot({
        'patterncnv': args.patterncnv,
        'r': args.Rscript,
        'gsize': args.gsize,
        'samtools': args.samtools,
        'bedtools': args.bedtools,
        'targets': args.exbaits,
        'refflat': args.exbaits,
        'ref': args.ref,
        'params': {
            'b': 10,
            'm': 20,
            'z': 1000,
            'x': 100,
            's': True
        }
    })

class BamSort(Proc):
    """Sort bam files

    Input:
        bamfile: The bam files to sort

    Output:
        outfile: The sorted bamfile

    Args:
        tool: Which tool to be used to do the sorting
        samtools: Path to `samtools`
        ncores: Number of cores to use
        by: Sort by `coord` or `name`
        params: Other params for the tool
    """
    lang = args.python
    input_keys = 'bamfile:file'
    output = 'outfile:file:{{in.bamfile | stem | append: ".sorted.bam"}}'
    script = f'file://{SCRIPT_DIR}/{MODULE}/BamSort.py'
    args = Diot({
        'tool': 'samtools',
        'samtools': args.samtools,
        'ncores': 1,
        'ref': args.ref,
        'by': 'coord', # or name
        'params': {}
    })