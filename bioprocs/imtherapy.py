"""A set of immunotherapy-related bioinformatics tools"""
from pyppl import Proc, Box
from . import params, delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pTopiary():
	"""
	@description:
		Predict mutation-derived cancer T-cell epitopes from somatic variants, tumor RNA expression data, and patient HLA type.
	@input:
		infile: The input VCF file.
			- Can be annotated with vcf_expression_annotator, as well as VEP
		afile: The HLA alleles in format of "HLA-A*02:01,HLA-B*07:02" or a file with one per line.
			- You can have star in them: "HLA-A*02:01,HLA-B*07:02"
	@output:
		outfile: The output file.
			- A html file will also be generated with suffix '.html'
		outdir: The output directory with output file and all intermediate files
	@args:
		topiary       (path): Path to topiary.
		netmhc        (path): Path to netmhc.
		netmhcpan     (path): Path to netmhcpan.
		netmhciipan   (path): Path to netmhciipan.
		netmhccons    (path): Path to netmhccons.
		smm           (path): Path to smm.
		smm_pmbec     (path): Path to smm-pmbec.
		mhc_predictor (path): The MHC binding predictor. Could be one of:
			- netmhc: Local NetMHC predictor (Topiary will attempt to automatically detect whether NetMHC 3.x or 4.0 is available)
			- netmhcpan: Local NetMHCpan predictor
			- netmhciipan: Local NetMHCIIpan predictor
			- netmhccons: Local NetMHCcons
			- random: Random IC50 values
			- smm: Local SMM predictor
			- smm-pmbec: Local SMM-PMBEC predictor
			- netmhcpan-iedb: Use NetMHCpan via the IEDB web API
			- netmhccons-iedb: Use NetMHCcons via the IEDB web API
			- smm-iedb: Use SMM via the IEDB web API
			- smm-pmbec-iedb: Use SMM-PMBEC via the IEDB web API
		params (Box): Other parameters for topiary
		tmpdir (str): Temporary directory for running local MHC predictors
	"""
	return Box(
		desc   = 'Epitope prediction using Topiary',
		lang   = params.python.value,
		input  = 'infile:file, afile:var',
		output = [
			'outfile:file:{{i.infile | stem2}}.topiary/{{i.infile | stem2}}.txt',
			'outdir:dir:{{i.infile | stem2}}.topiary'],
		args   = Box(
			topiary       = params.topiary.value,
			netmhc        = params.netmhc.value,
			netmhcpan     = params.netmhcpan.value,
			netmhciipan   = params.netmhciipan.value,
			netmhccons    = params.netmhccons.value,
			genome        = params.genome.value,
			smm           = params.smm.value,
			smm_pmbec     = params.smm_pmbec.value,
			tmpdir        = params.tmpdir.value,
			refall        = params.refall.value,
			mhc_predictor = 'netmhcpan',
			params        = Box()
		)
	)

@procfactory
def _pVACseq():
	"""
	@description:
		pVACseq is a cancer immunotherapy pipeline for the identification of personalized Variant Antigens by Cancer Sequencing (pVACseq) that integrates tumor mutation and expression data (DNA- and RNA-Seq). It enables cancer immunotherapy research by using massively parallel sequence data to predicting tumor-specific mutant peptides (neoantigens) that can elicit anti-tumor T cell immunity. It is being used in studies of checkpoint therapy response and to identify targets for personalized cancer vaccines and adoptive T cell therapies. For more general information.
		See: http://www.genomemedicine.com/content/8/1/11
	"""
	return Box(
		desc   = 'Identification of Variant Antigens using tumor mutation and expression data.',
		lang   = params.python.value,
		input  = 'infile:file, afile:var',
		output = [
			'outfile:file:{{i.infile | stem}}.pvacseq/{{i.infile | stem}}.epitopes.txt',
			'outdir:dir:{{i.infile | stem}}.pvacseq'],
		args   = Box(
			# used to extract sample name
			bcftools   = params.bcftools.value,
			pvacseq    = params.pvacseq.value,
			allele     = '',
			bdtool     = 'netmhc',
			netmhc     = params.netmhc.value,
			iedb_mhc_i = params.iedb_mhc_i.value,
			nthread    = 1,
			params     = Box(epitope_length = 9)
		)
	)

@procfactory
def _pOptiType():
	"""
	@description:
		Precision HLA typing from next-generation sequencing data using OptiType
		See: https://github.com/FRED-2/OptiType
	@input:
		fqfile1: The bam file or the fastq file of the 1st pair-end
		fqfile2: The fastq file of the 2nd pair-end
			- We will try to extract hla sequences from the bam/fastq file.
	@output:
		outfile: The output file with HLA types
		outdir: The output directory
	@args:
		optitype (str): Path to OptiTypePipeline.py
		picard   (str): Path to picard, used to convert bam files to fastq files.
			- Samtools not used, as it may mark paired reads with different names, which causes bwa crash.
		bwa      (str): Path to bwa
		hlaref   (str): Path to HLA reference file
		params   (Box): Other parameters for optitype
		nthread  (int): Number of threads to use
	@requires:
		[OptiType](https://github.com/FRED-2/OptiType)
		bwa
		picard
	"""
	return Box(
		desc   = 'Precision HLA typing from next-generation sequencing data using OptiType',
		lang   = params.python.value,
		input  = 'fqfile1:file, fqfile2:file',
		output = [
			'outfile:file:{% 	assign outdir = i.fqfile1, i.fqfile2 | \
								? :bool(_[1]) | \
								= :__import__("os").path.commonprefix | \
								! :_[0] | stem2 | .rstrip: "._[]" | \
								@append: ".optitype" %}{{outdir}}/{{outdir}}.txt',
			'outdir:dir:{{		i.fqfile1, i.fqfile2 | \
								? :bool(_[1]) | \
								= :__import__("os").path.commonprefix | \
								! :_[0] | stem2 | .rstrip: "._[]"}}.optitype'],
		args = Box(
			optitype = params.optitype.value,
			picard   = params.picard.value,
			bwa      = params.bwa.value,
			hlaref   = params.hlaref.value,
			nthread  = params.nthread.value,
			params   = Box(dna = True, bwa = Box())
		)
	)

@procfactory
def _pHLA_LA():
	"""
	@input:
		infile: The bam file
	@output:
		outfile: The output file with HLA types
		outdir: The output directory
	@args:
		hla_la   (path): Path to HLA-LA.pl
		picard   (path): Path to picard
		bwa      (path): Path to bwa
		samtools (path): Path to samtools
		nthread  (int) : Number of threads to use
		params   (Box) : Other parameters for `HLA-LA.pl`
	@requires:
		[HLA-LA](https://github.com/DiltheyLab/HLA-LA)
	"""
	return Box(
		desc   = 'HLA typing using HLA-LA',
		lang   = params.python.value,
		input  = 'infile:file',
		output = [	'outfile:file:{{i.infile | stem2}}/{{i.infile | stem2}}.hlala.txt',
					'outdir:dir:{{i.infile | stem2}}'],
		args   = Box(
			hla_la   = params.hla_la.value,
			picard   = params.picard.value,
			bwa      = params.bwa.value,
			samtools = params.samtools.value,
			java     = params.java.value,
			params   = Box(),
			nthread  = 1,
		)
	)

@procfactory
def _pVcf2ProteinSeq():
	"""
	@input:
		infile: The input VCF file
	@output:
		outdir: The output directory containing the protein sequences
	"""
	return Box(
		desc = 'Generate protein sequences from VCF file using pVACseq',
		lang   = params.python.value,
		input = 'infile:file',
		output = 'outdir:dir:{{i.infile | stem2}}.protseqs',
		args = Box(
			pvacseq = params.pvacseq.value,
			lens    = [15, 17, 19, 21],
			params  = Box(),
			nthread = 1,
		)
	)

@procfactory
def _pNetMHC():
	"""
	@input:
		infile: Peptide sequences in fasta format or list of peptide sequences, optionally with 2nd columns of numeric values.
		afile: A list (separated by comma) or a file of alleles, one per line. For example:
			- HLA-A*03:79,HLA-A*04:03
	@output:
		outfile: The output file of binding affinity of all combinations (Peptide ~ HLA allele)
	@args:
		netmhc  (path)        : Path to netMHC
		isfa    (bool)        : Whether the input file is fasta sequence file (otherwise it is peptides)
		params  (Box)         : Other parameters for netMHC
		lens    (int|str|list): Peptide length
		tmpdir  (path)        : Temporary directory for netMHC
		nthread (int)         : Number of threads to use.
			- netMHC itself does not multi-threading.
			- We will split the input file and put them into multiple threads
	"""
	return Box(
		desc   = 'Run the netmhc to predict binding affinity between HLA-allele and peptides',
		input  = 'infile:file, afile:var',
		output = 'outfile:file:{{i.infile | stem2}}.netmhc.xls',
		lang   = params.python.value,
		args   = Box(
			netmhc  = params.netmhc.value,
			isfa    = True,
			params  = Box(v = True),
			lens    = 9,
			tmpdir  = params.tmpdir.value,
			nthread = 1
		)
	)