from pyppl import Aggr
from bioprocs.common import pFile2Proc
from bioprocs.rnaseq import pExprDir2Matrix, pRNASeqDEG, pExprPlot
from bioprocs.marray import pCELDir2Matrix, pMArrayDEG
from bioprocs.resource import pTxt
from bioprocs.gsea import pExprMat2GCT, pSampleinfo2CLS, pGSEA, pEnrichr

"""
@name:
	aExpDir2Deg
@description:
	From expfiles to degs with sample info file.
@depends:
	```
	pExprDir2Matrix[*] \
	                      pRNASeqDEG[!]
	pSampleInfo[*]     /
	```
"""
aExpDir2Deg = Aggr(
	pExprDir2Matrix,
	pFile2Proc.copy(id = 'pSampleInfo'),
	pRNASeqDEG,
	depends = False
)
# Dependences
aExpDir2Deg.starts             = aExpDir2Deg.pExprDir2Matrix, aExpDir2Deg.pSampleInfo
aExpDir2Deg.ends               = aExpDir2Deg.pRNASeqDEG
aExpDir2Deg.pRNASeqDEG.depends = aExpDir2Deg.pExprDir2Matrix, aExpDir2Deg.pSampleInfo
# Input
aExpDir2Deg.pRNASeqDEG.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aExpDir2Deg.pExprDir2Matrix.args.boxplot  = True
aExpDir2Deg.pExprDir2Matrix.args.heatmap  = True
aExpDir2Deg.pExprDir2Matrix.args.histplot = True
aExpDir2Deg.pRNASeqDEG.args.maplot       = True
aExpDir2Deg.pRNASeqDEG.args.heatmap      = True


"""
@name:
	aExpDir2DegGSEA
@description:
	From expfiles to degs with sample info file and do GSEA.
@depends:
	```
	                  /   pRNASeqDEG[!]  -> pEnrichr[!]
	pExprDir2Matrix[*] \
	                      pExprMat2GCT     \
	pSampleInfo[*]     /                   \
	                   \  pSampleinfo2CLS --  pGSEA[!]
	pGMTFetcher[*]     ____________________/
	```
"""
aExpDir2DegGSEA = Aggr(
	pExprDir2Matrix,
	pFile2Proc.copy(id = 'pSampleInfo'),
	pTxt.copy(id = 'pGMTFetcher'),
	pExprPlot,
	pExprMat2GCT,
	pSampleinfo2CLS,
	pGSEA,
	pRNASeqDEG,
	pEnrichr,
	depends = False
)
# Default input:
aExpDir2DegGSEA.pGMTFetcher.input = ['KEGG_2016_gmt']
# Dependences
aExpDir2DegGSEA.starts                  = aExpDir2DegGSEA.pExprDir2Matrix,  \
										  aExpDir2DegGSEA.pSampleInfo,     \
										  aExpDir2DegGSEA.pGMTFetcher
aExpDir2DegGSEA.ends                    = aExpDir2DegGSEA.pRNASeqDEG,      \
										  aExpDir2DegGSEA.pGSEA,           \
										  aExpDir2DegGSEA.pEnrichr
aExpDir2DegGSEA.pExprMat2GCT.depends     = aExpDir2DegGSEA.pExprDir2Matrix
aExpDir2DegGSEA.pSampleinfo2CLS.depends = aExpDir2DegGSEA.pSampleInfo
aExpDir2DegGSEA.pGSEA.depends           = aExpDir2DegGSEA.pExprMat2GCT,     \
										  aExpDir2DegGSEA.pSampleinfo2CLS, \
										  aExpDir2DegGSEA.pGMTFetcher
aExpDir2DegGSEA.pRNASeqDEG.depends      = aExpDir2DegGSEA.pExprDir2Matrix,  \
										  aExpDir2DegGSEA.pSampleInfo
aExpDir2DegGSEA.pEnrichr.depends        = aExpDir2DegGSEA.pRNASeqDEG
# Input
aExpDir2DegGSEA.pGSEA.input      = lambda ch1, ch2, ch3: \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aExpDir2DegGSEA.pRNASeqDEG.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aExpDir2DegGSEA.pExprDir2Matrix.args.boxplot  = True
aExpDir2DegGSEA.pExprDir2Matrix.args.heatmap  = True
aExpDir2DegGSEA.pExprDir2Matrix.args.histplot = True
aExpDir2DegGSEA.pRNASeqDEG.args.maplot       = True
aExpDir2DegGSEA.pRNASeqDEG.args.heatmap      = True
aExpDir2DegGSEA.pGMTFetcher.args.header      = False

"""
@name:
	aRnaseqExpMat2DegGSEA
@description:
	From expfiles to degs with sample info file and do GSEA.
@depends:
	```
	                  /   pRNASeqDEG[!]  -> pEnrichr[!]
	pExpMat[*]        \
	                      pExprMat2GCT     \
	pSampleInfo[*]     /                   \
	                   \  pSampleinfo2CLS --  pGSEA[!]
	pGMTFetcher[*]     ____________________/
	```
"""
aRnaseqExpMat2DegGSEA = Aggr(
	pFile2Proc.copy(id = 'pExpMat'),
	pFile2Proc.copy(id = 'pSampleInfo'),
	pTxt.copy(id = 'pGMTFetcher'),
	pExprMat2GCT,
	pSampleinfo2CLS,
	pGSEA,
	pRNASeqDEG,
	pEnrichr,
	depends = False
)
# Default input:
aRnaseqExpMat2DegGSEA.pGMTFetcher.input = ['KEGG_2016_gmt']
# Dependences
aRnaseqExpMat2DegGSEA.starts                  = aRnaseqExpMat2DegGSEA.pExpMat,  \
												aRnaseqExpMat2DegGSEA.pSampleInfo,     \
												aRnaseqExpMat2DegGSEA.pGMTFetcher
aRnaseqExpMat2DegGSEA.ends                    = aRnaseqExpMat2DegGSEA.pRNASeqDEG,      \
												aRnaseqExpMat2DegGSEA.pGSEA,           \
												aRnaseqExpMat2DegGSEA.pEnrichr
aRnaseqExpMat2DegGSEA.pExprMat2GCT.depends     = aRnaseqExpMat2DegGSEA.pExpMat
aRnaseqExpMat2DegGSEA.pSampleinfo2CLS.depends = aRnaseqExpMat2DegGSEA.pSampleInfo
aRnaseqExpMat2DegGSEA.pGSEA.depends           = aRnaseqExpMat2DegGSEA.pExprMat2GCT,     \
												aRnaseqExpMat2DegGSEA.pSampleinfo2CLS, \
												aRnaseqExpMat2DegGSEA.pGMTFetcher
aRnaseqExpMat2DegGSEA.pRNASeqDEG.depends      = aRnaseqExpMat2DegGSEA.pExpMat,  \
												aRnaseqExpMat2DegGSEA.pSampleInfo
aRnaseqExpMat2DegGSEA.pEnrichr.depends        = aRnaseqExpMat2DegGSEA.pRNASeqDEG
# Input
aRnaseqExpMat2DegGSEA.pGSEA.input      = lambda ch1, ch2, ch3: \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aRnaseqExpMat2DegGSEA.pRNASeqDEG.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aRnaseqExpMat2DegGSEA.pRNASeqDEG.args.maplot       = True
aRnaseqExpMat2DegGSEA.pRNASeqDEG.args.heatmap      = True
aRnaseqExpMat2DegGSEA.pGMTFetcher.args.header      = False



"""
@name:
	aCELDir2DEG
@description:
	From celfils to degs with sample info file and do GSEA.
@procs:
	`pCELDir2Matrix[*]`: Convert CEL files to matrix
	`pSampleInfo[*]`   : Receive sample information file. Copied from `pFile2Proc`
	`pGMTFetcher[*]`   : Download the GMT file for GSEA. Copied from `pTxt`
	`pExprPlot[!]`     : Plot the expression distributions.
	`pExprMat2GCT`     : Convert expression matrix to GCT file for GSEA
	`pGSEA[!]`         : Do GSEA
	`pMArrayDEG[!]`    : Call DEGs
	`pEnrichr[!]`      : Do enrichment with DEGs
@modules
	`exprplot`: Plot expression distributions from expression matrix
		- `pCELDir2Matrix[*]` -> `pExprPlot[!]`
	`gsea`: Do GSEA
		- `pCELDir2Matrix[*]` -> `pExprMat2GCT`
		- `pSampleInfo[*]`    -> `pSampleinfo2CLS`
		- `pExprMat2GCT`, `pCELDir2Matrix[*]`, `pGMTFetcher` -> `pGSEA[!]`
	`enrichr`: Do enrichr
		- `pCELDir2Matrix[*]`, `pSampleInfo[*]` -> `pMArrayDEG[!]` -> `pEnrichr[!]`
"""
aCELDir2DEG = Aggr(
	pCELDir2Matrix,
	pFile2Proc.copy(id = 'pSampleInfo'),
	pTxt.copy(id = 'pGMTFetcher'),
	pExprPlot,
	pExprMat2GCT,
	pSampleinfo2CLS,
	pGSEA,
	pMArrayDEG,
	pEnrichr,
	depends = False
)
# Defaults
aCELDir2DEG.starts = 'pCELDir2Matrix, pSampleInfo'
aCELDir2DEG.ends   = 'pMArrayDEG'

aCELDir2DEG['pGMTFetcher'].input        = ['KEGG_2016_gmt']
aCELDir2DEG['pMArrayDEG'].args.maplot   = True
aCELDir2DEG['pMArrayDEG'].args.heatmap  = True
aCELDir2DEG['pGMTFetcher'].args .header = False
# Depends
aCELDir2DEG['pMArrayDEG', ].depends = ['pCELDir2Matrix, pSampleInfo']
# Input adjustments
aCELDir2DEG['pGSEA'].input = lambda ch1, ch2, ch3:       \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aCELDir2DEG['pMArrayDEG'].input = lambda ch1, ch2:       \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Configs
aCELDir2DEG.module('exprplot', depends = {
	'pExprPlot': 'pCELDir2Matrix'
}, ends = 'pExprPlot')
aCELDir2DEG.module('gsea', 'pGMTFetcher', {
	'pGSEA'          : 'pExprMat2GCT, pSampleinfo2CLS, pGMTFetcher',
	'pExprMat2GCT'   : 'pCELDir2Matrix',
	'pSampleinfo2CLS': 'pSampleInfo'
}, 'pGSEA')
aCELDir2DEG.module('enrichr', depends = {
	'pEnrichr': 'pMArrayDEG'
}, ends = 'pEnrichr')
# Enable all functions
aCELDir2DEG.on()

