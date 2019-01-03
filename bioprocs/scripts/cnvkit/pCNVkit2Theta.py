from os import path, makedirs
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

{% python from os import path %}
cnsfile  = {{i.cnsfile | quote}}
cnnfile  = {{i.cnnfile | quote}}
snvfile  = {{i.snvfile | quote}}
outdir   = {{o.outdir | quote}}
nthread  = {{args.nthread | quote}}
ckparams = {{args.ckparams | repr}}
ttparams = {{args.ttparams | repr}}
cnvkit   = {{args.cnvkit | quote}}
theta    = {{args.theta | quote}}
theta_in = {{job.outdir | path.join: "THetA2_input" | quote}}
theta_ic = {{job.outdir | path.join: "THetA2_input", fn2(i.cnsfile) + '.theta2.interval_count' | quote}}
theta_ts = {{job.outdir | path.join: "THetA2_input", fn(i.cnsfile) + '.normal.snp_formatted.txt' | quote}}
theta_ns = {{job.outdir | path.join: "THetA2_input", fn(i.cnsfile) + '.tumor.snp_formatted.txt' | quote}}
prefix   = {{o.outdir | path.join: fn2(i.cnsfile) | quote}}

if not path.exists(theta_in):
	makedirs(theta_in)

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

# region THetA example
# Generates Sample_T.BEST.results:
# /path/to/theta2/bin/RunTHetA Sample_T.interval_count \
#     --TUMOR_FILE Sample_T.tumor.snp_formatted.txt \
#     --NORMAL_FILE Sample_T.normal.snp_formatted.txt \
#     --BAF --NUM_PROCESSES `nproc` --FORCE
# endregion

# region THetA usage
# usage: RunTHetA.py [-h] [--TUMOR_FILE TUMOR_FILE] [--NORMAL_FILE NORMAL_FILE]
#                    [-n N] [-k K] [-t TAU] [-d DIR] [-p PRE] [-m MAX_NORMAL]
#                    [--NUM_PROCESSES NUM_PROCESSES]
#                    [--NUM_INTERVALS NUM_INTERVALS] [--BOUND_HEURISTIC BH]
#                    [--NORMAL_BOUND_HEURISTIC NBH] [--HEURISTIC_LB LB]
#                    [--HEURISTIC_UB UB] [--BOUNDS_ONLY] [--NO_MULTI_EVENT]
#                    [--RESULTS filename] [--FORCE] [--GET_VALUES]
#                    [--NO_INTERVAL_SELECTION] [--READ_DEPTH_FILE FILENAME]
#                    [--GRAPH_FORMAT GRAPH_FORMAT] [--BAF]
#                    [--RATIO_DEV RATIO_DEV] [--MIN_FRAC MIN_FRAC]
#                    [--NO_CLUSTERING]
#                    QUERY_FILE

# positional arguments:
#   QUERY_FILE            Interval file

# optional arguments:
#   -h, --help            show this help message and exit
#   --TUMOR_FILE TUMOR_FILE
#                         File containing allelic counts for tumor sample SNPs.
#   --NORMAL_FILE NORMAL_FILE
#                         File containing allelic counts for normal samlpe SNPs.
#   -n N, --N N           Number of subpopulations
#   -k K, --MAX_K K       The maximum value allowed for entries in C
#   -t TAU, --TAU TAU     Expected number of copies in normal genome
#   -d DIR, --DIR DIR     Directory where result file is written to
#   -p PRE, --OUTPUT_PREFIX PRE
#                         Prefix for output files created. By default, it will
#                         be the beginning of the input file name (i.e.if input
#                         filename were example.input, output filed would be
#                         example.output and example.withbounds)
#   -m MAX_NORMAL, --MAX_NORMAL MAX_NORMAL
#                         The maximum fraction to consider for normal. Only
#                         enforced for n=2
#   --NUM_PROCESSES NUM_PROCESSES
#                         The number of processes to be used
#   --NUM_INTERVALS NUM_INTERVALS
#                         The maximum number of intervals used by automatic
#                         interval selection.
#   --BOUND_HEURISTIC BH
#   --NORMAL_BOUND_HEURISTIC NBH
#   --HEURISTIC_LB LB
#   --HEURISTIC_UB UB
#   --BOUNDS_ONLY
#   --NO_MULTI_EVENT
#   --RESULTS filename
#   --FORCE
#   --GET_VALUES
#   --NO_INTERVAL_SELECTION
#   --READ_DEPTH_FILE FILENAME
#   --GRAPH_FORMAT GRAPH_FORMAT
#                         Options are .pdf, .jpg, .png, .eps
#   --BAF                 Option to run the BAF model.
#   --RATIO_DEV RATIO_DEV
#                         The deviation away from 1.0 that a ratio must be to be
#                         considered a potential copy number aberration.
#   --MIN_FRAC MIN_FRAC   The minimum fraction of the genome that must have a
#                         potential copy number aberration to be a valid sample
#                         for THetA analysis.
#   --NO_CLUSTERING       Option to run THetA without clustering.
# endregion

ckparams[""] = cnsfile
ckparams.o   = theta_ic
if cnnfile:
	ckparams.r = cnnfile
if snvfile:
	ckparams.v = snvfile
if 'i' in ckparams:
	ckparams.i = {{args.ckparams, render | *&: 'i' in a | *:b(a['i']) if c else None | repr}}
if 'n' in ckparams:
	ckparams.n = {{args.ckparams, render | *&: 'n' in a | *:b(a['n']) if c else None | repr}}
runcmd('cd {theta_in!r}; {cnvkit} export theta {args}'.format(theta_in = theta_in, cnvkit = cnvkit, args = cmdargs(ckparams, equal = ' ')))

ttparams[""] = theta_ic
ttparams.NUM_PROCESSES = nthread
ttparams.GRAPH_FORMAT  = '.png'
ttparams.OUTPUT_PREFIX = prefix
if path.isfile(theta_ts):
	ttparams.TUMOR_FILE = theta_ts
if path.isfile(theta_ns):
	ttparams.NORMAL_FILE = theta_ns
runcmd('{theta} {args} >/dev/stderr'.format(outdir = outdir, theta = theta, args = cmdargs(ttparams, equal = ' ')))
