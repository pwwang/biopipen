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
		infile: The input VCF or MAF file.
		mhcallele: The MHC alleles in format of "HLA-A02:01,HLA-B07:02" or a file with one per line.
			- You can have star in them: "HLA-A*02:01,HLA-B*07:02"
		gexpr: Cufflinks FPKM tracking file containing gene expression estimates.
		texpr: Cufflinks FPKM tracking file containing transcript expression estimates.
	@output:
		outfile: The output csv file.
			- A html file will also be generated with suffix '.html'
	@args:
		topiary       (str): Path to topiary.
		netmhc        (str): Path to netmhc.
		netmhcpan     (str): Path to netmhcpan.
		netmhciipan   (str): Path to netmhciipan.
		netmhccons    (str): Path to netmhccons.
		smm           (str): Path to smm.
		smm_pmbec     (str): Path to smm-pmbec.
		mhc_predictor (str): The MHC binding predictor. Could be one of:
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
		wildtype (bool): Also output wildtype peptide binding affinity.
			- Only available for local MHC binding affinity predictor
		tmpdir (str): Temporary directory for running local MHC predictors
	"""
	return Box(
		desc   = 'Epitope prediction using Topiary',
		lang   = params.python.value,
		input  = 'infile:file, mhcallele:var, gexpr:file, texpr:file',
		output = 'outfile:file:{{i.infile | stem}}.txt',
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
			mhc_predictor = 'netmhcpan',
			wildtype      = False,
			params        = Box({
				'output-csv-sep': '\t',

			})
		)
	)