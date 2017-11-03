from pyppl import Aggr
from bioprocs.common import pFile2Proc
from bioprocs.rnaseq import pExpdir2Matrix, pRnaseqDeg
from bioprocs.marray import pCeldir2Matrix, pMarrayDeg
from bioprocs.resource import pTxt
from bioprocs.gsea import pExpmat2Gct, pSampleinfo2Cls, pGSEA, pEnrichr

"""
@name:
	aExpDir2Deg
@description:
	From expfiles to degs with sample info file.
@depends:
	```
	pExpdir2Matrix[*] \
	                      pRnaseqDeg[!]
	pSampleInfo[*]     /
	```
"""
aExpDir2Deg = Aggr(
	pExpdir2Matrix,
	pFile2Proc.copy(newid = 'pSampleInfo'),
	pRnaseqDeg,
	depends = False
)
# Dependences
aExpDir2Deg.starts             = aExpDir2Deg.pExpdir2Matrix, aExpDir2Deg.pSampleInfo
aExpDir2Deg.ends               = aExpDir2Deg.pRnaseqDeg
aExpDir2Deg.pRnaseqDeg.depends = aExpDir2Deg.pExpdir2Matrix, aExpDir2Deg.pSampleInfo
# Input
aExpDir2Deg.pRnaseqDeg.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aExpDir2Deg.pExpdir2Matrix.args.boxplot  = True
aExpDir2Deg.pExpdir2Matrix.args.heatmap  = True
aExpDir2Deg.pExpdir2Matrix.args.histplot = True
aExpDir2Deg.pRnaseqDeg.args.maplot       = True
aExpDir2Deg.pRnaseqDeg.args.heatmap      = True


"""
@name:
	aExpDir2DegGSEA
@description:
	From expfiles to degs with sample info file and do GSEA.
@depends:
	```
	                  /   pRnaseqDeg[!]  -> pEnrichr[!]
	pExpdir2Matrix[*] \
	                      pExpmat2Gct     \
	pSampleInfo[*]     /                   \
	                   \  pSampleinfo2Cls --  pGSEA[!]
	pGMTFetcher        ____________________/
	```
"""
aExpDir2DegGSEA = Aggr(
	pExpdir2Matrix,
	pFile2Proc.copy(newid = 'pSampleInfo'),
	pTxt.copy(newid = 'pGMTFetcher'),
	pExpmat2Gct,
	pSampleinfo2Cls,
	pGSEA,
	pRnaseqDeg,
	pEnrichr,
	depends = False
)
# Default input:
aExpDir2DegGSEA.pGMTFetcher.input = ['KEGG_2016_gmt']
# Dependences
aExpDir2DegGSEA.starts                  = aExpDir2DegGSEA.pExpdir2Matrix,  \
										  aExpDir2DegGSEA.pSampleInfo,     \
										  aExpDir2DegGSEA.pGMTFetcher
aExpDir2DegGSEA.ends                    = aExpDir2DegGSEA.pRnaseqDeg,      \
										  aExpDir2DegGSEA.pGSEA,           \
										  aExpDir2DegGSEA.pEnrichr
aExpDir2DegGSEA.pExpmat2Gct.depends     = aExpDir2DegGSEA.pExpdir2Matrix
aExpDir2DegGSEA.pSampleinfo2Cls.depends = aExpDir2DegGSEA.pSampleInfo
aExpDir2DegGSEA.pGSEA.depends           = aExpDir2DegGSEA.pExpmat2Gct,     \
										  aExpDir2DegGSEA.pSampleinfo2Cls, \
										  aExpDir2DegGSEA.pGMTFetcher
aExpDir2DegGSEA.pRnaseqDeg.depends      = aExpDir2DegGSEA.pExpdir2Matrix,  \
										  aExpDir2DegGSEA.pSampleInfo
aExpDir2DegGSEA.pEnrichr.depends        = aExpDir2DegGSEA.pRnaseqDeg
# Input
aExpDir2DegGSEA.pGSEA.input      = lambda ch1, ch2, ch3: \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aExpDir2DegGSEA.pRnaseqDeg.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aExpDir2DegGSEA.pExpdir2Matrix.args.boxplot  = True
aExpDir2DegGSEA.pExpdir2Matrix.args.heatmap  = True
aExpDir2DegGSEA.pExpdir2Matrix.args.histplot = True
aExpDir2DegGSEA.pRnaseqDeg.args.maplot       = True
aExpDir2DegGSEA.pRnaseqDeg.args.heatmap      = True
aExpDir2DegGSEA.pGMTFetcher.args.header      = False

"""
@name:
	aRnaseqExpMat2DegGSEA
@description:
	From expfiles to degs with sample info file and do GSEA.
@depends:
	```
	                  /   pRnaseqDeg[!]  -> pEnrichr[!]
	pExpMat[*]        \
	                      pExpmat2Gct     \
	pSampleInfo[*]     /                   \
	                   \  pSampleinfo2Cls --  pGSEA[!]
	pGMTFetcher[*]     ____________________/
	```
"""
aRnaseqExpMat2DegGSEA = Aggr(
	pFile2Proc.copy(newid = 'pExpMat'),
	pFile2Proc.copy(newid = 'pSampleInfo'),
	pTxt.copy(newid = 'pGMTFetcher'),
	pExpmat2Gct,
	pSampleinfo2Cls,
	pGSEA,
	pRnaseqDeg,
	pEnrichr,
	depends = False
)
# Default input:
aRnaseqExpMat2DegGSEA.pGMTFetcher.input = ['KEGG_2016_gmt']
# Dependences
aRnaseqExpMat2DegGSEA.starts                  = aRnaseqExpMat2DegGSEA.pExpMat,  \
												aRnaseqExpMat2DegGSEA.pSampleInfo,     \
												aRnaseqExpMat2DegGSEA.pGMTFetcher
aRnaseqExpMat2DegGSEA.ends                    = aRnaseqExpMat2DegGSEA.pRnaseqDeg,      \
												aRnaseqExpMat2DegGSEA.pGSEA,           \
												aRnaseqExpMat2DegGSEA.pEnrichr
aRnaseqExpMat2DegGSEA.pExpmat2Gct.depends     = aRnaseqExpMat2DegGSEA.pExpMat
aRnaseqExpMat2DegGSEA.pSampleinfo2Cls.depends = aRnaseqExpMat2DegGSEA.pSampleInfo
aRnaseqExpMat2DegGSEA.pGSEA.depends           = aRnaseqExpMat2DegGSEA.pExpmat2Gct,     \
												aRnaseqExpMat2DegGSEA.pSampleinfo2Cls, \
												aRnaseqExpMat2DegGSEA.pGMTFetcher
aRnaseqExpMat2DegGSEA.pRnaseqDeg.depends      = aRnaseqExpMat2DegGSEA.pExpMat,  \
												aRnaseqExpMat2DegGSEA.pSampleInfo
aRnaseqExpMat2DegGSEA.pEnrichr.depends        = aRnaseqExpMat2DegGSEA.pRnaseqDeg
# Input
aRnaseqExpMat2DegGSEA.pGSEA.input      = lambda ch1, ch2, ch3: \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aRnaseqExpMat2DegGSEA.pRnaseqDeg.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aRnaseqExpMat2DegGSEA.pRnaseqDeg.args.maplot       = True
aRnaseqExpMat2DegGSEA.pRnaseqDeg.args.heatmap      = True
aRnaseqExpMat2DegGSEA.pGMTFetcher.args.header      = False

"""
@name:
	aCelDir2Deg
@description:
	From CEL files to degs with sample info file.
@depends:
	```
	pCeldir2Matrix[*] \
	                      pMarrayDeg[!]
	pSampleInfo[*]     /
	```
"""
aCelDir2Deg = Aggr(
	pCeldir2Matrix,
	pFile2Proc.copy(newid = 'pSampleInfo'),
	pMarrayDeg,
	depends = False
)
# Dependences
aCelDir2Deg.starts             = aCelDir2Deg.pCeldir2Matrix, aCelDir2Deg.pSampleInfo
aCelDir2Deg.ends               = aCelDir2Deg.pMarrayDeg
aCelDir2Deg.pMarrayDeg.depends = aCelDir2Deg.pCeldir2Matrix, aCelDir2Deg.pSampleInfo
# Input
aCelDir2Deg.pMarrayDeg.input = lambda ch1, ch2: \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aCelDir2Deg.pCeldir2Matrix.args.boxplot  = True
aCelDir2Deg.pCeldir2Matrix.args.heatmap  = True
aCelDir2Deg.pCeldir2Matrix.args.histplot = True
aCelDir2Deg.pMarrayDeg.args.maplot       = True
aCelDir2Deg.pMarrayDeg.args.heatmap      = True

"""
@name:
	aCelDir2DegGSEA
@description:
	From celfils to degs with sample info file and do GSEA.
@depends:
	```
	                  /   pMarrayDeg[!]  -> pEnrichr[!]
	pCeldir2Matrix[*] \
	                      pExpmat2Gct     \
	pSampleInfo[*]     /                   \
	                   \  pSampleinfo2Cls --  pGSEA[!]
	pGMTFetcher        ____________________/
	```
"""
aCelDir2DegGSEA = Aggr(
	pCeldir2Matrix,
	pFile2Proc.copy(newid = 'pSampleInfo'),
	pTxt.copy(newid = 'pGMTFetcher'),
	pExpmat2Gct,
	pSampleinfo2Cls,
	pGSEA,
	pMarrayDeg,
	pEnrichr,
	depends = False
)
# Default input:
aCelDir2DegGSEA.pGMTFetcher.input = ['KEGG_2016_gmt']
# Dependences
aCelDir2DegGSEA.starts                  = aCelDir2DegGSEA.pCeldir2Matrix, \
										  aCelDir2DegGSEA.pSampleInfo,    \
										  aCelDir2DegGSEA.pGMTFetcher
aCelDir2DegGSEA.ends                    = aCelDir2DegGSEA.pMarrayDeg,     \
										  aCelDir2DegGSEA.pGSEA,          \
										  aCelDir2DegGSEA.pEnrichr
aCelDir2DegGSEA.pExpmat2Gct.depends     = aCelDir2DegGSEA.pCeldir2Matrix
aCelDir2DegGSEA.pSampleinfo2Cls.depends = aCelDir2DegGSEA.pSampleInfo
aCelDir2DegGSEA.pGSEA.depends           = aCelDir2DegGSEA.pExpmat2Gct,    \
										  aCelDir2DegGSEA.pSampleinfo2Cls,\
										  aCelDir2DegGSEA.pGMTFetcher
aCelDir2DegGSEA.pMarrayDeg.depends      = aCelDir2DegGSEA.pCeldir2Matrix, \
										  aCelDir2DegGSEA.pSampleInfo
aCelDir2DegGSEA.pEnrichr.depends        = aCelDir2DegGSEA.pMarrayDeg
# Input
aCelDir2DegGSEA.pGSEA.input      = lambda ch1, ch2, ch3: \
	ch1.repRow(ch2.length() * ch3.length())              \
	   .cbind(ch2.repRow(ch1.length() * ch3.length()))   \
	   .cbind(ch3.repRow(ch1.length() * ch2.length())) 
aCelDir2DegGSEA.pMarrayDeg.input = lambda ch1, ch2:      \
	ch1.colAt(0).repRow(ch2.length()).cbind(ch2.repRow(ch1.length()))
# Args
aCelDir2DegGSEA.pCeldir2Matrix.args.boxplot  = True
aCelDir2DegGSEA.pCeldir2Matrix.args.heatmap  = True
aCelDir2DegGSEA.pCeldir2Matrix.args.histplot = True
aCelDir2DegGSEA.pMarrayDeg.args.maplot       = True
aCelDir2DegGSEA.pMarrayDeg.args.heatmap      = True
aCelDir2DegGSEA.pGMTFetcher.args.header      = False
