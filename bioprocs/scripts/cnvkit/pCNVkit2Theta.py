from os import path
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

{% python from os import path %}
cnsfile  = {{i.cnsfile | quote}}
cnnfile  = {{i.cnnfile | quote}}
outfile  = {{o.outfile | quote}}
nthread  = {{args.nthread | quote}}
params   = {{args.params | repr}}
cnvkit   = {{args.cnvkit | quote}}

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

openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)

params[""] = cnsfile
params.o   = path.basename(outfile)
if cnnfile:
	params.r = cnnfile

runcmd('cd {outdir!r}; {nthr}; {cnvkit} export theta {args}'.format(
	nthr   = openblas_nthr,
	outdir = path.dirname(outfile),
	cnvkit = cnvkit,
	args   = cmdargs(params, equal = ' ')
))



