"""Operations on TCGA MAF file"""
from pyppl import Proc
from diot import Diot
from .utils import fs2name
# Next-step processing after VCF file being generated.
from . import params, proc_factory
from .tumhet import pTMBurden

pGTMatAddRs = proc_factory(
	desc = "Add rs id to a genotype matrix",
	config = Diot(annotate = """
	@name:
		pGTMatAddRs
	@description:
		Add rs id to a genotype matrix
	@input:
		infile: The input genotype matrix, columns are samples, rows are mutations in format:
			- `<chr>_<pos>_<ref>_<alt>` or `<chr>_<pos>_<name>_<ref>_<alt>`
			- has to be sorted by coordinates.
			- `chromsomes` have to be in order of `args.chrorder`
	@output:
		outfile: The output genotype matrix. Row names will turn into:
			- `<chr>_<pos>_<rs>_<ref>_<alt>`
	@args:
		`dbsnp`: the dbsnp vcf file used to annotation the snps.
			- assume sorted by coordinates
			- `chromsomes` have to be in order of `args.chrorder`
		`notfound`: What to used if RS id not found. Default: `NOVEL`
			- `None/Fase` to skip to record
		`exist`: What if RS id exists? Default: `keep`
			- `keep`: Keep the RS ID and skip seeking
			- `force`: Force using the RS ID being found to replace the old one.
		`chrorder`: The chromsome order. Default: `<params.chrorder>`
	"""))
pGTMatAddRs.input         = 'infile:file'
pGTMatAddRs.output        = 'outfile:file:{{i.infile | bn}}'
pGTMatAddRs.args.dbsnp    = params.dbsnp_all.value
pGTMatAddRs.args.chrorder = params.chrorder.value
pGTMatAddRs.args.notfound = 'NOVEL'
pGTMatAddRs.args.exist    = 'keep'
pGTMatAddRs.lang          = params.python.value

pGTMat2Plink = proc_factory(
	desc = 'Convert a genotype matrix to plink binary files',
	config = Diot(annotate = """
	@name:
		pGTMat2Plink
	@description:
		Convert a genotype matrix to plink binary files
	@input:
		infile: The genotype matrix, probably generated by `pVcf2GTMat`
		`metafile:file`: The metadata file.
			- column names could be `['FID', 'IID', 'PID', 'MID', 'Sex', 'Pheno']`, see plink's `ped` format
			- row names are samples
	@output:
		`outdir:dir`: The output directory. Default: `{{i.infile | fn}}.plink`
	@args:
		`plink`:   The path to `plink`
		`keeptxt`: Keep the text files (.ped and .map) or not. Default: `False`
	@requires:
		`plink 1.x`
	"""))
pGTMat2Plink.input        = 'infile:file, metafile:file'
pGTMat2Plink.output       = 'outdir:dir:{{i.infile | fn}}.plink'
pGTMat2Plink.args.plink   = params.plink.value
pGTMat2Plink.args.keeptxt = False
pGTMat2Plink.args.chrmaps = {'X': 23, 'Y': 24, 'XY': 25, 'M': 26, 'MT': 26}
pGTMat2Plink.lang         = params.python.value

pGTMat2Bed = proc_factory(
	desc = 'Convert a genotype matrix to bed file',
	config = Diot(annotate = """
	@name:
		pGTMat2Bed
	@description:
		Convert a genotype matrix to a bed file containing the coordinates of the mutations
	@input:
		infile: The genotype matrix. Row names must follow `<chr>_<pos>_<rsid>_<ref>_<alt>`
	@output:
		outfile: The output bed file. Default: `outfile:file:{{i.infile | fn}}.bed`
	@args:
		`ncol`  : How many columns of bed to output. Default: `6`.
			- Possible values: 3, 6, 8, 66, 88
			- If `8`, then reference and alternative alleles will be 7th and 8th column
			- If `66`, then genotypes will be attached to BED6
			- If `88`, then genotypes will be attached to `ncol = 8`
		`name`  : Use the neat name (usually rsid) or full name (row names). Default: `neat`
		`inopts`: Options to read the input file. Default: `Diot(cnames = True)`
	"""))
pGTMat2Bed.input       = 'infile:file'
pGTMat2Bed.output      = 'outfile:file:{{i.infile | fn}}.bed'
pGTMat2Bed.args.inopts = Diot(cnames = True)
pGTMat2Bed.args.ncol   = 6
pGTMat2Bed.args.name   = 'neat' # full
pGTMat2Bed.lang        = params.python.value

pCallRate = proc_factory(
	desc   = 'Calculate sample/snp call rate from single sample vcfs',
	lang   = params.Rscript.value,
	config = Diot(annotate = """
	@name:
		pCallRate
	@description:
		Calculate sample/snp call rate from single sample vcfs
	@input:
		indir:     The dir containing the vcfs
	@output:
		`outsample:file`: The report of call rate for each sample
		`figsample:file`: The bar chat of sample call rates
		`outsnp:file`:    The report of call rate for each snp
		`figsnp:file`:    The bar chat of snp call rates
	"""))
pCallRate.input            = "indir:file"
pCallRate.output           = "outdir:dir:{{ i.indir | fn }}.callrate"
pCallRate.args.histplotggs = []
pCallRate.args.devpars     = Diot({'res':300, 'width':2000, 'height':2000})

pCepip = proc_factory(
	desc = 'Run cepip for input mutations.',
	config = Diot(annotate = """
	@name:
		pCepip
	@description:
		Run CEPIP.
	@input:
		infile: The input file (vcf or avinput)
	@output:
		outfile: The cepip result file
	@args:
		`cepip`:    The path of cepip
		`cell` :    The related cell line
		`params`:   Other params for cepip
	@requires:
		[`cepip`](http://jjwanglab.org/cepip/)
	"""))
pCepip.input       = "infile:file"
pCepip.output      = "outfile:file:{{i.infile | fn}}.cepip.txt"
pCepip.args.cepip  = params.cepip.value
pCepip.args.cell   = ""
pCepip.args.params = Diot()
pCepip.lang        = params.python.value

pMutSig = proc_factory(
	desc = 'Run MutSig.',
	config = Diot(annotate = """
	@name:
		pMutSig
	@description:
		MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.
		For more information, see Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013).
		See [dcumentation](http://archive.broadinstitute.org/cancer/cga/mutsig_run)
	@input:
		infile: mutation table
	@output:
		`outdir:dir`: The output directory
	@args:
		`mutsig` : The path to `run_MutSigCV.sh`, default: 'mutsig'
		`mcr`    : The Matlab MCR path
		`cvrg`   : coverage table
		`cvrt`   : covariates table
		`mutdict`: mutation_type_dictionary_file
		`chrdir` : chr_files_hg18 or chr_files_hg19
	@requires:
		[MutSig](http://archive.broadinstitute.org/cancer/cga/mutsig_download)
	"""))
pMutSig.input        = 'infile:file'
pMutSig.output       = "outdir:dir:{{i.infile | fn}}.mutsig"
pMutSig.args.cvrg    = params.mutsig_cvrg.value
pMutSig.args.cvrt    = params.mutsig_cvrt.value
pMutSig.args.mutdict = params.mutsig_mutdict.value
pMutSig.args.chrdir  = params.mutsig_chrdir.value
pMutSig.args.mutsig  = params.mutsig.value
pMutSig.args.mcr     = params.mcr.value

pMafLiftover = proc_factory(
	desc = 'Liftover a maf file from one assembly to another',
	config = Diot(annotate = """
	@name:
		pMafLiftover
	@description:
		Liftover maf file from one assembly to another
	@input:
		infile: The input maf file
	@output:
		outfile: The output maf file
	@args:
		`liftover`: The liftOver program.
		`lochain`:  The liftOver chain file.
		`genome`:   The target genome.
	@requires:
		liftOver from UCSC
	"""))
pMafLiftover.input         = 'infile:file'
pMafLiftover.output        = 'outfile:file:{{i.infile | fn | lambda x: x if x.endswith(".maf") else x + ".maf"}}'
pMafLiftover.args.liftover = params.liftover.value
pMafLiftover.args.lochain  = params.lochain.value
pMafLiftover.args.genome   = params.genome.value
pMafLiftover.lang          = params.python.value

pMafMerge = proc_factory(
	desc = 'Merge maf files',
	config = Diot(annotate = """
	@input:
		`infiles:files`: The maf files
	@output:
		outfile: The merged maf file
	@args:
		`excols`: How to deal with extra columns other than 34 standard columns from TCGA.
			- merge(default): Merge the columns, if one not exists, fill with __UNKNOWN__.
			- discard: Just discard the extra columns, with only 34 columns left. So you can also put just one maf file in the indir with some columns missed to fill it with standard columns.
	"""),
	input  = 'infiles:files',
	output = 'outfile:file:{{i.infiles | fs2name}}.maf',
	lang   = params.python.value,
	args   = Diot(excols  = 'merge'),  # discard,
	envs   = Diot(fs2name = fs2name)
)

pMaf2Mat = proc_factory(
	desc = 'Convert maf file to a gene-based mutation matrix',
	config = Diot(annotate = """
	@name:
		pMaf2Mat
	@description:
		Convert maf file to a gene(row)-sample(column) matrix
	@input:
		infile: The input file
	@output:
		outfile: The output matrix
	@args:
		`mutypes`: Provide manual list of variant classifications to be counted, only effective when `args.binary = False`. Default: `None` (all counted)
		`binary` : Just generate a binary matrix instead of a count matrix. Default: `False`
		`na`: What value to use for no mutations reported on a gene. Default: `0`
		`samfn`  : A function (in r) to transform the sample names. Default: `function(sample) sample`
	"""))
pMaf2Mat.input        = 'infile:file'
pMaf2Mat.output       = 'outfile:file:{{i.infile | fn}}.mat.txt'
pMaf2Mat.args.binary  = False
pMaf2Mat.args.mutypes = None
pMaf2Mat.args.na      = 0
pMaf2Mat.args.samfn   = 'function(sample) sample'
pMaf2Mat.lang         = params.Rscript.value

pMaftools = proc_factory(
	desc   = 'Basic analysis on somatic mutations of a group of samples',
	config = Diot(annotate = """
	@input:
		indir: The input directory or a single maf file. A directory could contain:
			- `*.maf` or `*.maf.gz` file (required)
			- `*.annot.tsv` or `*.annot.txt` file (see: https://github.com/PoisonAlien/maftools/blob/master/inst/extdata/tcga_laml_annot.tsv)
			- `all_lesions.conf_*.txt`: Gistic cnv data
			- `amp_genes.conf_*.txt`: Gistic cnv data
			- `del_genes.conf_*.txt`: Gistic cnv data
			- `scores.gistic`: Gistic cnv data
			- `*.seg.txt`: CBS segments data
		msdir: MutSig result directory.
			- Must have`*sig_genes.txt` to do pancancer comparison.
	@output:
		outdir: The output directory
	@args:
		ngenes : Top number of genes to plot for some plots.
		extypes: Exclude mutation types, only consider nonsynonymous mutations.
		isTCGA : If the maf file is from TCGA?
		ref    : The reference file for signature plot.
		plot   : Which plots to plot.
		params : The extra parameters for each plot function.
		devpars: The parameters for plot device.
		nthread: Number of threads used for multiple plot of one type.
	@requires:
		[Maftools](https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html)
	"""),
	input  = 'indir:file, msdir:file',
	output = 'outdir:dir:{{i.indir | fn}}.maftools',
	lang   = params.Rscript.value,
	args   = Diot(
		ngenes  = 10,
		isTCGA  = False,
		genome  = params.genome.value,
		devpars = Diot(res = 300, height = 2000, width = 2000),
		nthread = 1,
		# exclude mutation types, only consider nonsynonymous mutations.
		extypes = ["Intron", "5'UTR", "3'UTR", "IGR", "5'Flank", "3'Flank", "Silent"],
		plot    = Diot(
			summary        = True,
			oncoplot       = True,
			oncostrip      = True,
			titv           = True,
			lollipop       = True,
			cbsseg         = True,
			rainfall       = True,
			tcgacomp       = True,
			vaf            = True,
			genecloud      = True,
			gisticGenome   = True,
			gisticBubble   = True,
			gisticOncoplot = True,
			somInteraction = True,
			oncodrive      = True,
			pfam           = True,
			pancan         = True,
			survival       = True,
			heterogeneity  = True,
			signature      = True,
		),
		params = Diot(
			summary        = Diot(rmOutlier = True, addStat = 'median', dashboard = True),
			oncoplot       = Diot(),
			oncostrip      = Diot(),
			titv           = Diot(),
			lollipop       = Diot(),
			cbsseg         = Diot(labelAll = True),
			rainfall       = Diot(detectChangePoints = True),
			tcgacomp       = Diot(),
			vaf            = Diot(flip = True),
			genecloud      = Diot(minMut = 3),
			gisticGenome   = Diot(markBands = 'all'),
			gisticBubble   = Diot(),
			gisticOncoplot = Diot(),
			somInteraction = Diot(),
			oncodrive      = Diot(minMut = 5, pvalMethod = 'zscore', fdrCutOff = 0.1, useFraction = True),
			pfam           = Diot(),
			pancan         = Diot(qval = 0.1, label = 1),
			survival       = Diot(),
			heterogeneity  = Diot(),
			signature      = Diot(nTry = 6, plotBestFitRes = False),
		)
	)
)

pMutationSigs = proc_factory(
	desc = 'Find similar COSMIC mutation signatures for MAF file.',
	config = Diot(annotate = """
	@name:
		pMutationSigs
	@description:
		Find similar COSMIC mutation signatures for MAF file
		using https://github.com/pwwang/deconstruct_sigs_py
	@input:
		infile: The input maf file.
	@output:
		`outdir:dir`: The output directory
	@args:
		`font_family`: Font family for plotting.
		`font_weight`: Font weight for plotting.
		`sig_cutoff` : Significance cutoff for signatures.
		`err_thres`  : The threshold to top the iteration.
		`ref`        : The reference genome.
	"""))
pMutationSigs.input            = 'infile:file'
pMutationSigs.output           = 'outdir:dir:{{i.infile | fn2}}.signature'
pMutationSigs.args.font_family = 'Arial'
pMutationSigs.args.font_weight = 'bold'
pMutationSigs.args.sig_cutoff  = 0.05
pMutationSigs.args.err_thres   = 1e-3
pMutationSigs.args.ref         = params.ref.value
pMutationSigs.lang             = params.python.value

pSnpEff = proc_factory(
	desc   = 'Discovery and characterization of artifactual mutations',
	config = Diot(annotate = """
	@name:
		pSnpEff
	@description:
		This is the default command. It is used for annotating variant filed (e.g. VCF files).
	@input:
		infile:  The input file
	@output:
		`outdir:file`: The directory containing output anntated file, snpEff_genes.txt and snpEff_summary.html
	@args:
		`snpEff`:       The snpEff executable, default: "snpEff"
		`params`:    Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"
		`genome`:    The genome used for annotation, default: "hg19"
		`informat`:  The format of input file [vcf or bed], default: "vcf"
		`outformat`: The format of output file [vcf, gatk, bed, bedAnn], default: "vcf"
		`csvStats`:  Whether to generate csv stats file, default: True.
		`htmlStats`: Whether to generate the html summary file, default: False.
		`javamem`:   The memory to use. Default: '-Xms1g -Xmx8g'
	@requires:
		[snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)
	"""))
pSnpEff.input  = "infile:file"
pSnpEff.output = "outdir:dir:{{infile | fn}}.snpeff"
pSnpEff.args   = { "snpEff": "snpEff", "javamem": "-Xms1g -Xmx8g", "genome": "hg19", "informat": "vcf", "outformat": "vcf", "csvStats": True, "htmlStats": False, "params": "" }
pSnpEff.script = """
csvfile="{{outdir}}/{{infile | fn}}.csvstat"
sumfile="{{outdir}}/{{infile | fn}}.html"
outfile="{{outdir}}/{{infile | fn}}.snpEff.vcf"
csvStats=""
if [[ "{{args.csvStats}}" == "True" ]]; then
	csvStats="-csvStats \\"$csvfile\\""
fi
stats=""
if [[ "{{args.htmlStats}}" == "True" ]]; then
	stats="-stats \\"$sumfile\\""
fi
echo {{args.snpEff}} {{args.javamem}} -i {{args.informat}} -o {{args.outformat}} $csvStats $stats {{args.params}} {{args.genome}} "{{infile}}"
{{args.snpEff}} {{args.javamem}} -i {{args.informat}} -o {{args.outformat}} $csvStats $stats {{args.params}} {{args.genome}} "{{infile}}" > "$outfile"
"""

pDToxoG = proc_factory(
	desc  = 'Run D-ToxoG on MAF files',
	config = Diot(annotate = """
	@description:
		Run D-ToxoG on MAF files with columns listed at:
		https://software.broadinstitute.org/cancer/cga/dtoxog
		However,`i_picard_oxoQ` is required but not documented. This process
		will try to put 0's for the field to let the program run.
		To be stricter, use picard's CollectOxoGMetrics to get this column.
	@input:
		infile: The input maf file
	@output:
		outfile: The output maf file.
			- with extra columns added:
			- `pox` --  p-value that the call is actually an artifact
			- `qox` --  false detection rate score.
			- `pox_cutoff` --  minimum pox score for artifact.
			- `isArtifactMode`:
				- Variant is C>A, G>T: 1
				- Variant is not C>A or G>T: 0
			- `oxoGCut`:
				- Variant is marked as artifact: 1
				- Variant is not an artifact: 0
	@args:
		dtoxog : D-ToxoG executable.
		nthread: Maximum threads used by matlab.
		keep   : Whether keep those artifact mode mutations in output MAF file or not.
		params : Other parameters for `startFilterMAFFile`
			- See more in `startFilterMAFFile.m` or run `dtoxog` directly
	"""),
	input  = 'infile:file',
	output = 'outfile:file:{{i.infile | stem}}.dtoxog.maf',
	lang   = params.python.value,
	args   = Diot(
		nthread = 1,
		keep    = True,
		dtoxog  = params.dtoxog.value,
		params  = Diot(isGeneratingPlots = True, globalPoxoG = .96, artifactThresholdRate = .01)
	))

pMaf2Vcf = proc_factory(
	desc   = 'Convert MAF file back to VCF files',
	config = Diot(annotate = """
	@input:
		infile: The input MAF file
	@output:
		outfile: Output multi-sample VCF containing all TN-pairs
		outdir: A directory with output VCF files
	@args:
		maf2vcf (str): Path to `maf2vcf.pl` from `vcf2maf` tool.
		merge: Whether merge the VCF for individual samples or not.
		params: Other parameters for `maf2vcf.pl`
	"""),
	lang   = params.python.value,
	input  = 'infile:file',
	output = [
		'outfile:file:{{i.infile | stem}}.vcfs/{{i.infile | stem}}.vcf',
		'outdir:dir:{{i.infile | stem}}.vcfs'],
	args = Diot(
		maf2vcf = params.maf2vcf.value,
		ref     = params.ref.value,
		params  = Diot({'per-tn-vcfs': True})),
)


pMafAddChr = proc_factory(
	desc = 'Add chr to chromosome if not present',
	config = Diot(annotate = """
	@input:
		infile: The input MAF file
	@output:
		outfile: The output MAF file
	"""),
	lang = params.python.value,
	input = 'infile:file',
	output = 'outfile:file:{{i.infile | stem}}.maf'
)

pMafExtractSample = proc_factory(
	desc   = 'Filter MAF file with given samples',
	config = Diot(annotate = """
	@input:
		infile: The input MAF file
		samfile: The tumor sample file
			- Could also be a list of samples, separated by comma
			- Overwrite `args.samples`
	@output:
		outfile: The output file contain only the given samples
	@args:
		samples (str|list): The samples or a list of samples
	"""),
	input  = 'infile:file, samfile:var',
	output = 'outfile:file:{{i.infile | stem}}.subset.maf',
	lang   = params.python.value,
	args   = Diot(samples = [])
)

pMafExtractClass = proc_factory(
	desc   = 'Filter MAF file with given classes',
	config = Diot(annotate = """
	@input:
		infile: The input MAF file
		classfile: The variant class file
			- Could also be a list of classes, separated by comma
			- Overwrite `args.classes`
	@output:
		outfile: The output file contain only the given classes
	@args:
		classes (str|list): The classes or a list of classes
	"""),
	input  = 'infile:file, classfile:var',
	output = 'outfile:file:{{i.infile | stem}}.subclass.maf',
	lang   = params.python.value,
	args   = Diot(classes = [])
)

pMafGetSamples = proc_factory(
	desc   = 'Get samples from MAF file',
	config = Diot(annotate = """
	@input:
		infile: The input MAF file
	@output:
		outfile: The output sample file
	@args:
		out: Which samples to output, tumor, normal or both.
		outopts: The output options
	"""),
	lang   = params.python.value,
	input  = 'infile:file',
	output = 'outfile:file:{{i.infile | stem}}.samples.txt',
	args   = Diot(out = 'both', outopts = Diot(cnames = True))
)

pMafFromTsv = proc_factory(
	desc   = 'Convert a TSV file with essential columns to a MAF file',
	config = Diot(annotate = """
	@description:
		Convert a TSV file with essential columns to a MAF file.
		See https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#protected-maf-file-structure for MAF file columns.
		Essential columns include: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2
		End_Position = Start_Position + 1 if not provided
	@input:
		infile: The TSV file, cnames(header) is required
	@output:
		output: The output MAF file
	@args:
		inopts   (Diot): Options to read input file.
		genome   (str) : The genome to fill the `NCBI_Build` column
		missing  (str) : The string used to fill other missing values.
		tumor    (str) : The string used to fill `Tumor_Sample_Barcode` if not provided.
		normal   (str) : The string used to fill `Matched_Norm_Sample_Barcode` if not provided.
		full     (bool): Output the basic MAF file (32 columns) or full (126 columns)
		refall   (path): The gene GTF file to extract gene names according to the coordinates
		bedtools (bool): Used to extract gene names
	"""),
	input  = 'infile:file',
	output = 'outfile:file:{{i.infile | stem}}.maf',
	lang   = params.python.value,
	args   = Diot(
		inopts   = Diot(cnames = True),
		genome   = params.genome.value,
		missing  = '__UNKNOWN__',
		tumor    = 'TUMOR',
		normal   = 'NORMAL',
		full     = True,
		refall   = params.refall.value,
		bedtools = params.bedtools.value,
	)
)
