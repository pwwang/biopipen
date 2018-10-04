from pyppl import Box
from deconstructSigs import DeconstructSigs, Settings

infile = {{i.infile | quote}}
outdir = {{o.outdir | quote}}
args   = {{args | repr}}

Settings.font_family = args.get('font_family', 'Arial')
Settings.font_weight = args.get('font_weight', 'bold')
Settings.sig_cutoff  = args.get('sig_cutoff' , .05)
Settings.err_thres   = args.get('err_thres'  , 1e-3)
ref = args.get('ref')

ds      = DeconstructSigs(maf = infile, hg19_fasta_path = ref, output_folder = outdir)
weights = ds.which_signatures()
ds.figures(weights, explanations=True)