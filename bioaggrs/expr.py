from pyppl import Aggr
from bioprocs.common import pFile2Proc
from bioprocs.rnaseq import pExprDir2Matrix, pRNAseqDEG
from bioprocs.marray import pCELDir2Matrix, pMArrayDEG
from bioprocs.resource import pTxt
from bioprocs.gsea import pExpmat2Gct, pSampleinfo2Cls, pGSEA, pEnrichr

"""
@name:
	aExpDir2Deg
@description:
	From expfiles to degs with sample info file.
@depends:
	```
	pExprDir2Matrix[*] \
	                      pRNAseqDEG[!]
	pSampleInfo[*]     /
	```
"""
aExpDir2Deg = Aggr(
	pExprDir2Matrix,
	pFile2Proc.copy(id = 'pSampleInfo'),
	pRNAseqDEG,
	depends = False
)
# Dependences
aExpDir2Deg.starts             = aExpDir2Deg.pExprDir2Matrix, aExpDir2Deg.pSampleInfo
aExpDir2Deg.ends               = aExpDir2Deg.pRNAseqDEG
aExpDir2Deg.pRNAseqDEG.depends = aExpDir2Deg.pExprDir2Matrix, aExpDir2Deg.pSampleInfo
# Input
aExpDir2Deg.pRNAseqDEG.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aExpDir2Deg.pExprDir2Matrix.args.boxplot  = True
aExpDir2Deg.pExprDir2Matrix.args.heatmap  = True
aExpDir2Deg.pExprDir2Matrix.args.histplot = True
aExpDir2Deg.pRNAseqDEG.args.maplot       = True
aExpDir2Deg.pRNAseqDEG.args.heatmap      = True


"""
@name:
	aExpDir2DegGSEA
@description:
	From expfiles to degs with sample info file and do GSEA.
@depends:
	```
	                  /   pRNAseqDEG[!]  -> pEnrichr[!]
	pExprDir2Matrix[*] \
	                      pExpmat2Gct     \
	pSampleInfo[*]     /                   \
	                   \  pSampleinfo2Cls --  pGSEA[!]
	pGMTFetcher[*]     ____________________/
	```
"""
aExpDir2DegGSEA = Aggr(
	pExprDir2Matrix,
	pFile2Proc.copy(id = 'pSampleInfo'),
	pTxt.copy(id = 'pGMTFetcher'),
	pExpmat2Gct,
	pSampleinfo2Cls,
	pGSEA,
	pRNAseqDEG,
	pEnrichr,
	depends = False
)
# Default input:
aExpDir2DegGSEA.pGMTFetcher.input = ['KEGG_2016_gmt']
# Dependences
aExpDir2DegGSEA.starts                  = aExpDir2DegGSEA.pExprDir2Matrix,  \
										  aExpDir2DegGSEA.pSampleInfo,     \
										  aExpDir2DegGSEA.pGMTFetcher
aExpDir2DegGSEA.ends                    = aExpDir2DegGSEA.pRNAseqDEG,      \
										  aExpDir2DegGSEA.pGSEA,           \
										  aExpDir2DegGSEA.pEnrichr
aExpDir2DegGSEA.pExpmat2Gct.depends     = aExpDir2DegGSEA.pExprDir2Matrix
aExpDir2DegGSEA.pSampleinfo2Cls.depends = aExpDir2DegGSEA.pSampleInfo
aExpDir2DegGSEA.pGSEA.depends           = aExpDir2DegGSEA.pExpmat2Gct,     \
										  aExpDir2DegGSEA.pSampleinfo2Cls, \
										  aExpDir2DegGSEA.pGMTFetcher
aExpDir2DegGSEA.pRNAseqDEG.depends      = aExpDir2DegGSEA.pExprDir2Matrix,  \
										  aExpDir2DegGSEA.pSampleInfo
aExpDir2DegGSEA.pEnrichr.depends        = aExpDir2DegGSEA.pRNAseqDEG
# Input
aExpDir2DegGSEA.pGSEA.input      = lambda ch1, ch2, ch3: \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aExpDir2DegGSEA.pRNAseqDEG.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aExpDir2DegGSEA.pExprDir2Matrix.args.boxplot  = True
aExpDir2DegGSEA.pExprDir2Matrix.args.heatmap  = True
aExpDir2DegGSEA.pExprDir2Matrix.args.histplot = True
aExpDir2DegGSEA.pRNAseqDEG.args.maplot       = True
aExpDir2DegGSEA.pRNAseqDEG.args.heatmap      = True
aExpDir2DegGSEA.pGMTFetcher.args.header      = False

"""
@name:
	aRnaseqExpMat2DegGSEA
@description:
	From expfiles to degs with sample info file and do GSEA.
@depends:
	```
	                  /   pRNAseqDEG[!]  -> pEnrichr[!]
	pExpMat[*]        \
	                      pExpmat2Gct     \
	pSampleInfo[*]     /                   \
	                   \  pSampleinfo2Cls --  pGSEA[!]
	pGMTFetcher[*]     ____________________/
	```
"""
aRnaseqExpMat2DegGSEA = Aggr(
	pFile2Proc.copy(id = 'pExpMat'),
	pFile2Proc.copy(id = 'pSampleInfo'),
	pTxt.copy(id = 'pGMTFetcher'),
	pExpmat2Gct,
	pSampleinfo2Cls,
	pGSEA,
	pRNAseqDEG,
	pEnrichr,
	depends = False
)
# Default input:
aRnaseqExpMat2DegGSEA.pGMTFetcher.input = ['KEGG_2016_gmt']
# Dependences
aRnaseqExpMat2DegGSEA.starts                  = aRnaseqExpMat2DegGSEA.pExpMat,  \
												aRnaseqExpMat2DegGSEA.pSampleInfo,     \
												aRnaseqExpMat2DegGSEA.pGMTFetcher
aRnaseqExpMat2DegGSEA.ends                    = aRnaseqExpMat2DegGSEA.pRNAseqDEG,      \
												aRnaseqExpMat2DegGSEA.pGSEA,           \
												aRnaseqExpMat2DegGSEA.pEnrichr
aRnaseqExpMat2DegGSEA.pExpmat2Gct.depends     = aRnaseqExpMat2DegGSEA.pExpMat
aRnaseqExpMat2DegGSEA.pSampleinfo2Cls.depends = aRnaseqExpMat2DegGSEA.pSampleInfo
aRnaseqExpMat2DegGSEA.pGSEA.depends           = aRnaseqExpMat2DegGSEA.pExpmat2Gct,     \
												aRnaseqExpMat2DegGSEA.pSampleinfo2Cls, \
												aRnaseqExpMat2DegGSEA.pGMTFetcher
aRnaseqExpMat2DegGSEA.pRNAseqDEG.depends      = aRnaseqExpMat2DegGSEA.pExpMat,  \
												aRnaseqExpMat2DegGSEA.pSampleInfo
aRnaseqExpMat2DegGSEA.pEnrichr.depends        = aRnaseqExpMat2DegGSEA.pRNAseqDEG
# Input
aRnaseqExpMat2DegGSEA.pGSEA.input      = lambda ch1, ch2, ch3: \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aRnaseqExpMat2DegGSEA.pRNAseqDEG.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aRnaseqExpMat2DegGSEA.pRNAseqDEG.args.maplot       = True
aRnaseqExpMat2DegGSEA.pRNAseqDEG.args.heatmap      = True
aRnaseqExpMat2DegGSEA.pGMTFetcher.args.header      = False

"""
@name:
	aCELDir2DEG
@description:
	From CEL files to degs with sample info file.
@depends:
	```
	pCELDir2Matrix[*] \
	                      pMArrayDEG[!]
	pSampleInfo[*]     /
	```
"""
aCELDir2DEG = Aggr(
	pCELDir2Matrix,
	pFile2Proc.copy(id = 'pSampleInfo'),
	pMArrayDEG,
	depends = False
)
# Dependences
aCELDir2DEG.starts             = aCELDir2DEG.pCELDir2Matrix, aCELDir2DEG.pSampleInfo
aCELDir2DEG.ends               = aCELDir2DEG.pMArrayDEG
aCELDir2DEG.pMArrayDEG.depends = aCELDir2DEG.pCELDir2Matrix, aCELDir2DEG.pSampleInfo
# Delegates
aCELDir2DEG.delegate('args', 'pMArrayDEG')
# Input
aCELDir2DEG.pMArrayDEG.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aCELDir2DEG.pCELDir2Matrix.args.boxplot  = True
aCELDir2DEG.pCELDir2Matrix.args.heatmap  = True
aCELDir2DEG.pCELDir2Matrix.args.histplot = True
aCELDir2DEG.pMArrayDEG.args.maplot       = True
aCELDir2DEG.pMArrayDEG.args.heatmap      = True

"""
@name:
	aCELDir2DEGGSEA
@description:
	From celfils to degs with sample info file and do GSEA.
@depends:
	```
	                  /   pMArrayDEG[!]  -> pEnrichr[!]
	pCELDir2Matrix[*] \
	                      pExpmat2Gct     \
	pSampleInfo[*]     /                   \
	                   \  pSampleinfo2Cls --  pGSEA[!]
	pGMTFetcher        ____________________/
	```
"""
aCELDir2DEGGSEA = Aggr(
	pCELDir2Matrix,
	pFile2Proc.copy(id = 'pSampleInfo'),
	pTxt.copy(id = 'pGMTFetcher'),
	pExpmat2Gct,
	pSampleinfo2Cls,
	pGSEA,
	pMArrayDEG,
	pEnrichr,
	depends = False
)
# Default input:
aCELDir2DEGGSEA.pGMTFetcher.input = ['KEGG_2016_gmt']
# Dependences
aCELDir2DEGGSEA.starts                  = aCELDir2DEGGSEA.pCELDir2Matrix, \
										  aCELDir2DEGGSEA.pSampleInfo,    \
										  aCELDir2DEGGSEA.pGMTFetcher
aCELDir2DEGGSEA.ends                    = aCELDir2DEGGSEA.pMArrayDEG,     \
										  aCELDir2DEGGSEA.pGSEA,          \
										  aCELDir2DEGGSEA.pEnrichr
aCELDir2DEGGSEA.pExpmat2Gct.depends     = aCELDir2DEGGSEA.pCELDir2Matrix
aCELDir2DEGGSEA.pSampleinfo2Cls.depends = aCELDir2DEGGSEA.pSampleInfo
aCELDir2DEGGSEA.pGSEA.depends           = aCELDir2DEGGSEA.pExpmat2Gct,    \
										  aCELDir2DEGGSEA.pSampleinfo2Cls,\
										  aCELDir2DEGGSEA.pGMTFetcher
aCELDir2DEGGSEA.pMArrayDEG.depends      = aCELDir2DEGGSEA.pCELDir2Matrix, \
										  aCELDir2DEGGSEA.pSampleInfo
aCELDir2DEGGSEA.pEnrichr.depends        = aCELDir2DEGGSEA.pMArrayDEG
# Input
aCELDir2DEGGSEA.pGSEA.input      = lambda ch1, ch2, ch3: \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aCELDir2DEGGSEA.pMArrayDEG.input = lambda ch1, ch2:      \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aCELDir2DEGGSEA.pCELDir2Matrix.args.boxplot  = True
aCELDir2DEGGSEA.pCELDir2Matrix.args.heatmap  = True
aCELDir2DEGGSEA.pCELDir2Matrix.args.histplot = True
aCELDir2DEGGSEA.pMArrayDEG.args.maplot       = True
aCELDir2DEGGSEA.pMArrayDEG.args.heatmap      = True
aCELDir2DEGGSEA.pGMTFetcher.args.header      = False
