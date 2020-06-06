from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.parallel import Parallel

infile  = Path({{i.infile | quote}})
outdir  = Path({{o.outdir | quote}})
pvacseq = {{args.pvacseq | quote}}
lens    = {{args.lens | list | repr}}
params  = {{ args.params | repr}}
nthread = {{args.nthread | repr}}

shell.load_config(pvacseq = pvacseq)

def do_one(protein_length):
	ps = params.copy()
	ps._ = [infile, protein_length, outdir.joinpath('%s.%s.peptide' % (
		infile.name.split('.')[0], protein_length))]
	shell.pvacseq.generate_protein_fasta(**ps).fg

para = Parallel(nthread = nthread)
para.run(do_one, [(plen,) for plen in lens])
