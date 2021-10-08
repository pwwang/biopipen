from .core.proc import Proc
from .core.args import args
from .core.defaults import SCRIPT_DIR

MODULE = 'fastx'

class TrimGalore(Proc):
    """Trim fastq files using trim_galore

    Input:
        fq1: fastq file for one pair-mate
        fq2: fastq file for the other

    Output:
        fq1out: The trimmed file for fq1
        fq2out: The trimmed file for fq2

    Args:
        trim_galore: Path to trim_galore
        params: Params for trim_galore
    """
    lang = args.python
    input_keys = 'fq1:file, fq2:file'
    output = [
        'fq1:file:{{in.fq1 | stem}}_trimmed.fq{{'
            'in.fq1 | .endswith: ".gz" ? ".gz" ! ""}}',
        'fq2:file:{{in.fq2 | stem}}_trimmed.fq{{'
            'in.fq1 | .endswith: ".gz" ? ".gz" ! ""}}'
    ]
    script = f'file://{SCRIPT_DIR}/{MODULE}/TrimGalore.py'
    args = {
        'trim_galore': args.trim_galore,
        'params': {
            'adapter': 'AGATCGGAAGAGC',
            'adapter2': 'AAATCAAAAAAAC',
            'paired': True
        }
    }
