from biopipen.utils.misc import run_command, dict_to_cli_args

infile = {{in.infile | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
outdir = {{out.outdir | quote}}  # pyright: ignore
perl = {{envs.perl | quote}}  # pyright: ignore
ref = {{envs.ref | repr}}  # pyright: ignore
samtools = {{envs.samtools | quote}}  # pyright: ignore
args = {{envs.args | repr}}  # pyright: ignore
maf2vcf = {{biopipen_dir | append: "/scripts/tcgamaf/maf2vcf.pl" | repr}}  # pyright: ignore

args['input-maf']  = infile
args['output-vcf'] = outfile
args['output-dir'] = outdir
args['ref-fasta']  = ref
args[''] = [perl, maf2vcf]

run_command(dict_to_cli_args(args, dashify=True), fg=True)
