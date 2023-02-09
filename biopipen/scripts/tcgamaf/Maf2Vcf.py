import cmdy

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

cmd = cmdy.perl(maf2vcf, _exe=perl, **args).hold()

print("Running:")
print(cmd.strcmd)

cmd.fg().run()
