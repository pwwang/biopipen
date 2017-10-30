from pyppl import Aggr
from bioprocs.common import pStr2File, pFile2Proc
from bioprocs.seq import pPromoters, pConsv, pConsvPerm
from bioprocs.bed import pBedGetfasta
from bioprocs.tfbs import pMotifScan

"""
@name:
	aTFBSOnPromoters
@description:
	Scan motifs on genes' promoter regions.
@depends:
	pPromoters[*] \
					pBedGetfasta \
								pMotifScan[!]
					pTFFile[*] /
"""
aTFBSOnPromoters = Aggr(
	pFile2Proc.copy(newid = 'pTFFile'),
	pPromoters,
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTFBSOnPromoters.pTFFile.runner    = 'local'
aTFBSOnPromoters.pPromoters.runner = 'local'
# delegate
aTFBSOnPromoters.delegate('args.up', 'pPromoters')
aTFBSOnPromoters.delegate('args.down', 'pPromoters')
aTFBSOnPromoters.delegate('args.genome', 'pPromoters')
aTFBSOnPromoters.delegate('args.ref', 'pBedGetfasta')
aTFBSOnPromoters.delegate('args.pval', 'pMotifScan')
aTFBSOnPromoters.delegate('args.tfmotifs', 'pMotifScan')
# depends
aTFBSOnPromoters.starts               = aTFBSOnPromoters.pTFFile,     aTFBSOnPromoters.pPromoters 
aTFBSOnPromoters.ends                 = aTFBSOnPromoters.pMotifScan
aTFBSOnPromoters.pMotifScan.depends   = aTFBSOnPromoters.pTFFile,     aTFBSOnPromoters.pBedGetfasta
aTFBSOnPromoters.pBedGetfasta.depends = aTFBSOnPromoters.pPromoters
# input
aTFBSOnPromoters.pMotifScan.input = lambda ch1, ch2: ch1.repRow(max(ch1.length(), ch2.length())).cbind(ch2)
# args
aTFBSOnPromoters.pBedGetfasta.args.params.name = True

"""
@name:
	aTFBSOnRegions
@description:
	Scan motifs on a given regions.
@depends:
	pBedGetfasta[*]  \
						pMotifScan[!]
		   pTFFile[*] /
"""
aTFBSOnRegions = Aggr(
	pFile2Proc.copy(newid = 'pTFFile'),
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTFBSOnRegions.pTFFile.runner     = 'local'
# delegate
aTFBSOnRegions.delegate('args.ref', 'pBedGetfasta')
aTFBSOnRegions.delegate('args.pval', 'pMotifScan')
# depends
aTFBSOnRegions.starts             = aTFBSOnRegions.pTFFile, aTFBSOnRegions.pBedGetfasta
aTFBSOnRegions.ends               = aTFBSOnRegions.pMotifScan
aTFBSOnRegions.pMotifScan.depends = aTFBSOnRegions.pTFFile, aTFBSOnRegions.pBedGetfasta
# args
aTFBSOnRegions.pBedGetfasta.args.params.name = True

"""
@name:
	aTFBSOnPromotersConsv
@description:
	Scan motifs on genes' promoter regions with conservation.
@depends:
	pPromoters[*] \
					pBedGetfasta \
								  pMotifScan
					   pTFFile[*] /	         \
					                          pConsv[!]
											 /
						    	pConsvPerm[*]
"""
aTFBSOnPromotersConsv = Aggr(
	pFile2Proc.copy(newid = 'pTFFile'),
	pPromoters,
	pConsvPerm,
	pBedGetfasta,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTFBSOnPromotersConsv.pTFFile.runner     = 'local'
aTFBSOnPromotersConsv.pPromoters.runner = 'local'
aTFBSOnPromotersConsv.pConsvPerm.input  = [0]
# delegate
aTFBSOnPromotersConsv.delegate('args.up', 'pPromoters')
aTFBSOnPromotersConsv.delegate('args.down', 'pPromoters')
aTFBSOnPromotersConsv.delegate('args.genome', 'pPromoters')
aTFBSOnPromotersConsv.delegate('args.ref', 'pBedGetfasta')
aTFBSOnPromotersConsv.delegate('args.pval', 'pMotifScan')
aTFBSOnPromotersConsv.delegate('args.threspval', 'pConsv')
aTFBSOnPromotersConsv.delegate('args.len', 'pConsvPerm')
aTFBSOnPromotersConsv.delegate('args.nperm', 'pConsvPerm')
aTFBSOnPromotersConsv.delegate('args.chrsizes', 'pConsvPerm')
# depends
aTFBSOnPromotersConsv.starts               = aTFBSOnPromotersConsv.pTFFile, aTFBSOnPromotersConsv.pPromoters, aTFBSOnPromotersConsv.pConsvPerm
aTFBSOnPromotersConsv.ends                 = aTFBSOnPromotersConsv.pConsv
aTFBSOnPromotersConsv.pConsv.depends       = aTFBSOnPromotersConsv.pMotifScan, aTFBSOnPromotersConsv.pConsvPerm
aTFBSOnPromotersConsv.pMotifScan.depends   = aTFBSOnPromotersConsv.pTFFile, aTFBSOnPromotersConsv.pBedGetfasta
aTFBSOnPromotersConsv.pBedGetfasta.depends = aTFBSOnPromotersConsv.pPromoters
# input
aTFBSOnPromotersConsv.pConsv.input = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTFBSOnPromotersConsv.pBedGetfasta.args.params.name = True
aTFBSOnPromotersConsv.pConsv.args.pval              = True


"""
@name:
	aTFBSOnRegionsConsv
@description:
	Scan motifs on a given regions with conservation score.
@depends:
	pBedGetfasta[*]  \
						pMotifScan
		   pTFFile[*] /             \
		                            pConsv[!]
				             	   /
					  pConsvPerm[*]
"""
aTFBSOnRegionsConsv = Aggr(
	pFile2Proc.copy(newid = 'pTFFile'),
	pBedGetfasta,
	pConsvPerm,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTFBSOnRegionsConsv.pTFFile.runner     = 'local'
aTFBSOnRegionsConsv.pConsvPerm.input  = [0]
# delegate
aTFBSOnRegionsConsv.delegate('args.ref', 'pBedGetfasta')
aTFBSOnRegionsConsv.delegate('args.pval', 'pMotifScan')
aTFBSOnRegionsConsv.delegate('args.threspval', 'pConsv')
aTFBSOnRegionsConsv.delegate('args.len', 'pConsvPerm')
aTFBSOnRegionsConsv.delegate('args.nperm', 'pConsvPerm')
aTFBSOnRegionsConsv.delegate('args.chrsizes', 'pConsvPerm')
# depends
aTFBSOnRegionsConsv.starts             = aTFBSOnRegionsConsv.pTFFile, aTFBSOnRegionsConsv.pBedGetfasta, aTFBSOnRegionsConsv.pConsvPerm
aTFBSOnRegionsConsv.ends               = aTFBSOnRegionsConsv.pConsv
aTFBSOnRegionsConsv.pConsv.depends     = aTFBSOnRegionsConsv.pMotifScan, aTFBSOnRegionsConsv.pConsvPerm
aTFBSOnRegionsConsv.pMotifScan.depends = aTFBSOnRegionsConsv.pTFFile, aTFBSOnRegionsConsv.pBedGetfasta
# input
aTFBSOnRegionsConsv.pConsv.input = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTFBSOnRegionsConsv.pBedGetfasta.args.params.name = True
aTFBSOnRegionsConsv.pConsv.args.pval              = True



