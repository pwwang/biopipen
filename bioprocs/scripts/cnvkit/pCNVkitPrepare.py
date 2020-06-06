"""Script for pCnvkitPrepare"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping,unsupported-assignment-operation
# pylint: disable=unsubscriptable-object

from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader
from bioprocs.utils.reference import bam_index, fa_index

infiles = {{i.infiles | repr}}
baits   = {{args.baits | repr}}
accfile = {{args.accfile | repr}}
ref     = {{args.ref | repr}}
params  = {{args.params | repr}}
prefix  = {{i.infiles | fs2name | quote}}
outdir  = {{job.outdir | quote}}
cnvkit  = {{args.cnvkit | quote}}
nthread = {{args.nthread | repr}}

fa_index(ref)

# get chroms
chroms = TsvReader(ref + '.fai', cnames=False).dump(0)

for infile in infiles:
    bam_index(infile)

shell.load_config(cnvkit = dict(
    _exe = cnvkit,
    _env = dict(
        OPENBLAS_NUM_THREADS = str(nthread),
        OMP_NUM_THREADS      = str(nthread),
        NUMEXPR_NUM_THREADS  = str(nthread),
        MKL_NUM_THREADS      = str(nthread)
    ),
    _cwd = outdir
))

# generate target file
if baits[-4:] in ('.gff', '.gtf'):
    #https://github.com/pwwang/pygff
    from gff import Gff
    baitfile = path.join(outdir, path.basename(baits[:-4]) + '.bait.bed')
    gff = Gff(baits)
    with open(baitfile, 'w') as f:
        for g in gff:
            # Exclude contigs not in reference file.
            if g['seqid'] not in chroms:
                continue
            f.write("{chrom}\t{start}\t{end}\t{name}\n".format(
                chrom = g['seqid'],
                start = g['start'],
                end   = g['end'],
                name  = g['attributes']['gene_id'],
            ))
    baits = baitfile

params_t   = params.target
params_t.o = path.join(outdir, prefix + '.bed')
shell.cnvkit.target(baits, **params_t).fg

# generate access file
if not accfile:
    accfile    = path.join(outdir, prefix + '.access.bed')
    params_a   = params.access
    params_a.o = accfile
    shell.cnvkit.access(ref, **params_a).fg

# autobin
params_b   = params.autobin
params_b.t = params_t.o
params_b.g = accfile
shell.cnvkit.autobin(*infiles, **params_b).fg
