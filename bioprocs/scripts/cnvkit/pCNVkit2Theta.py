from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell

cnsfile  = {{i.cnsfile | quote}}
cnnfile  = {{i.cnnfile | quote}}
outfile  = {{o.outfile | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params | repr}}
cnvkit   = {{args.cnvkit | quote}}

shell.load_config(cnvkit = dict(
	_exe = cnvkit,
	_env = dict(
		OPENBLAS_NUM_THREADS = str(nthread),
		OMP_NUM_THREADS      = str(nthread),
		NUMEXPR_NUM_THREADS  = str(nthread),
		MKL_NUM_THREADS      = str(nthread)
	),
	_cwd = path.dirname(outfile)
))

# region cnvkit export example
# cnvkit.py export theta Sample_T.cns reference.cnn -v Sample_Paired.vcf
# cnvkit.py export theta Sample_Tumor.cns Sample_Normal.cnr -o Sample.theta2.interval_count
# cnvkit.py export theta Sample_Tumor.cns -o Sample.theta2.interval_count
# endregion

# region cnvkit export theta usage
# usage: cnvkit.py export theta [-h] [-r REFERENCE] [-o OUTPUT] [-v VCF]
#                               [-i SAMPLE_ID] [-n NORMAL_ID]
#                               [-m MIN_VARIANT_DEPTH] [-z [ALT_FREQ]]
#                               tumor_segment

# positional arguments:
#   tumor_segment         Tumor-sample segmentation file from CNVkit (.cns).

# optional arguments:
#   -h, --help            show this help message and exit
#   -r REFERENCE, --reference REFERENCE
#                         Reference copy number profile (.cnn), or normal-sample
#                         bin-level log2 copy ratios (.cnr). Use if the
#                         tumor_segment input file does not contain a "weight"
#                         column.
#   -o OUTPUT, --output OUTPUT
#                         Output file name.

# To also output tables of SNP b-allele frequencies for THetA2:
#   -v VCF, --vcf VCF     VCF file containing SNVs observed in both the tumor
#                         and normal samples. Tumor sample ID should match the
#                         `tumor_segment` filename or be specified with -i
#                         /--sample-id.
#   -i SAMPLE_ID, --sample-id SAMPLE_ID
#                         Specify the name of the tumor sample in the VCF (given
#                         with -v/--vcf). [Default: taken the tumor_segment file
#                         name]
#   -n NORMAL_ID, --normal-id NORMAL_ID
#                         Corresponding normal sample ID in the input VCF.
#   -m MIN_VARIANT_DEPTH, --min-variant-depth MIN_VARIANT_DEPTH
#                         Minimum read depth for a SNP in the VCF to be counted.
#                         [Default: 20]
#   -z [ALT_FREQ], --zygosity-freq [ALT_FREQ]
#                         Ignore VCF's genotypes (GT field) and instead infer
#                         zygosity from allele frequencies. [Default if used
#                         without a number: 0.25]
# endregion

params.o   = path.basename(outfile)
if cnnfile:
	params.r = cnnfile

shell.fg.cnvkit.export('theta', cnsfile, **params)

