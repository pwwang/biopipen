"""Script for fastx.pFastqTrim"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping

from pathlib import Path
from diot import Diot
from bioprocs.utils import mem2, shell2 as shell

fq1          = {{i.fq1 | quote}}
fq2          = {{i.fq2 | quote}}
outfq1       = {{o.outfq1 | quote}}
outfq2       = {{o.outfq2 | quote}}
tool         = {{args.tool | .replace: ' ', '' | quote}}
cutadapt     = {{args.cutadapt | quote}}
skewer       = {{args.skewer |quote}}
trimmomatic  = {{args.trimmomatic | quote}}
params       = {{args.params | repr}}
nthread      = {{args.nthread | repr}}
gz           = {{args.gz | repr}}
mem          = {{args.mem | quote }}
minlen       = {{args.minlen |repr}}
minq         = {{args.minq |repr}}
cut5         = {{args.cut5 |repr}}
cut3         = {{args.cut3 |repr}}
adapter1     = {{args.adapter1 |quote}}
adapter2     = {{args.adapter2 |quote}}
bbmap_repair = {{args.bbmap_repair | quote}}
jobout       = Path({{job.outdir |quote}})

shell.load_config(
    cutadapt     = cutadapt,
    skewer       = skewer,
    trimmomatic  = trimmomatic,
    bbmap_repair = bbmap_repair,
)

def _seqrev (seq):
    d = {
        'A':'T', 'T':'A', 'G':'C', 'C':'G',
        'a':'t', 't':'a', 'g':'c', 'c':'g'
    }
    return ''.join([d[s] for s in seq])

def run_trimmomatic():
    """Run trimmomatic"""
    mem_java = mem2(mem, 'java').split()
    minlen2 = minlen * 2
    adfile = jobout / "adapters.fa"
    adfile.write_text("\n".join([
        ">PE1", _seqrev(adapter1),
        ">PE1_rc", adapter1,
        ">PE2", _seqrev(adapter2),
        ">PE2_rc", adapter2,
    ]))
    params2 = params.copy()
    params2.threads = nthread
    params2._prefix = '-'
    shell.trimmomatic(*mem_java, _sub=True).PE(
        params2, fq1, fq2, outfq1, "/dev/null",
        outfq2, "/dev/null",
        "ILLUMINACLIP:%s:2:30:10" % adfile,
        "LEADING:%s" % cut5,
        "TRAILING:%s" % cut3,
        "SLIDINGWINDOW:4:%s" % minq,
        "MINLEN:%s" % minlen2
    ).fg

def run_cutadapt():
    """Run cutadapt"""
    params2 = params.copy()
    params2.a = adapter1
    params2.A = adapter2
    params2.u = [cut5, 0-cut3]
    params2.U = [cut5, 0-cut3]
    params2.m = minlen
    params2.q = "{0}.{0}".format(minq)
    params2.o = outfq1
    params2.O = outfq2
    params2._ = [fq1, fq2]
    shell.cutadapt(**params2).fg

def run_skewer():
    """Run skewer"""
    params2 = params.copy()
    params2.m = 'pe'
    params2.t = nthread
    params2.x = adapter1
    params2.y = adapter2
    params2.Q = minq
    params2.l = minlen
    params2.z = gz
    params2.o = jobout / 'tmp'
    params2._ = [fq1, fq2]
    shell.skewer(**params2).fg
    out1 = jobout / 'tmp-trimmed-pair1.fastq'
    out2 = jobout / 'tmp-trimmed-pair2.fastq'
    if gz:
        out1 = str(out1) + '.gz'
        out2 = str(out2) + '.gz'
    shell.mv(out1, outfq1)
    shell.mv(out2, outfq2)

def run_bbmap_repair():
    """Run repair.sh from bbmap to correct some unmatched paired reads"""
    mem_java = mem2(mem, 'java').split()[1]
    unrepaired_fq1 = jobout.joinpath('unrepaired_' + Path(outfq1).name)
    unrepaired_fq2 = jobout.joinpath('unrepaired_' + Path(outfq2).name)
    shell.mv(outfq1, unrepaired_fq1)
    shell.mv(outfq2, unrepaired_fq2)

    params2 = params.copy()
    params2['in'] = unrepaired_fq1
    params2.in2 = unrepaired_fq2
    params2.out = outfq1
    params2.out2 = outfq2
    params2.overwrite = 'f'
    shell.bbmap_repair(mem_java, **params2).fg
    shell.rm_rf(unrepaired_fq1)
    shell.rm_rf(unrepaired_fq2)

TOOLS = {
    'trimmomatic': run_trimmomatic,
    'cutadapt'   : run_cutadapt,
    'skewer'     : run_skewer,
    'trimmomatic+bbmap_repair': lambda: run_trimmomatic() or run_bbmap_repair(),
    'cutadapt+bbmap_repair': lambda: run_cutadapt() or run_bbmap_repair(),
    'skewer+bbmap_repair': lambda: run_skewer() or run_bbmap_repair(),
}

try:
    TOOLS[tool]()
except KeyError:
    raise ValueError('{} not supported yet.'.format(tool))
