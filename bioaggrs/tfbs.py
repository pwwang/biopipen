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
					pMotif[*] /
"""
aTFBSOnPromoters = Aggr(
	pFile2Proc.copy(newid = 'pMotif'),
	pPromoters,
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTFBSOnPromoters.pMotif.runner     = 'local'
aTFBSOnPromoters.pPromoters.runner = 'local'
# delegate
aTFBSOnPromoters.delegate('args.mids', 'pMotifScan')
aTFBSOnPromoters.delegate('args.gnorm', 'pMotifScan')
aTFBSOnPromoters.delegate('args.up', 'pPromoters')
aTFBSOnPromoters.delegate('args.down', 'pPromoters')
aTFBSOnPromoters.delegate('args.genome', 'pPromoters')
aTFBSOnPromoters.delegate('args.ref', 'pBedGetfasta')
aTFBSOnPromoters.delegate('args.pval', 'pMotifScan')
# depends
aTFBSOnPromoters.starts               = aTFBSOnPromoters.pMotif,     aTFBSOnPromoters.pPromoters 
aTFBSOnPromoters.ends                 = aTFBSOnPromoters.pMotifScan
aTFBSOnPromoters.pMotifScan.depends   = aTFBSOnPromoters.pMotif,     aTFBSOnPromoters.pBedGetfasta
aTFBSOnPromoters.pBedGetfasta.depends = aTFBSOnPromoters.pPromoters
# args
aTFBSOnPromoters.pBedGetfasta.args.params.name = True
aTFBSOnPromoters.pMotifScan.args.gnorm         = True

"""
@name:
	aTFBSOnRegions
@description:
	Scan motifs on a given regions.
@depends:
	pBedGetfasta[*]  \
						pMotifScan[!]
		   pMotif[*] /
"""
aTFBSOnRegions = Aggr(
	pFile2Proc.copy(newid = 'pMotif'),
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTFBSOnRegions.pMotif.runner     = 'local'
# delegate
aTFBSOnRegions.delegate('args.mids', 'pMotifScan')
aTFBSOnRegions.delegate('args.gnorm', 'pMotifScan')
aTFBSOnRegions.delegate('args.ref', 'pBedGetfasta')
aTFBSOnRegions.delegate('args.pval', 'pMotifScan')
# depends
aTFBSOnRegions.starts             = aTFBSOnRegions.pMotif, aTFBSOnRegions.pBedGetfasta
aTFBSOnRegions.ends               = aTFBSOnRegions.pMotifScan
aTFBSOnRegions.pMotifScan.depends = aTFBSOnRegions.pMotif, aTFBSOnRegions.pBedGetfasta
# args
aTFBSOnRegions.pBedGetfasta.args.params.name = True
aTFBSOnRegions.pMotifScan.args.gnorm         = True

"""
@name:
	aTFBSOnPromotersConsv
@description:
	Scan motifs on genes' promoter regions with conservation.
@depends:
	pPromoters[*] \
					pBedGetfasta \
								  pMotifScan
					   pMotif[*] /	         \
					                          pConsv[!]
											 /
						    	pConsvPerm[*]
"""
aTFBSOnPromotersConsv = Aggr(
	pFile2Proc.copy(newid = 'pMotif'),
	pPromoters,
	pConsvPerm,
	pBedGetfasta,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTFBSOnPromotersConsv.pMotif.runner     = 'local'
aTFBSOnPromotersConsv.pPromoters.runner = 'local'
aTFBSOnPromotersConsv.pConsvPerm.input  = [0]
# delegate
aTFBSOnPromotersConsv.delegate('args.mids', 'pMotifScan')
aTFBSOnPromotersConsv.delegate('args.gnorm', 'pMotifScan')
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
aTFBSOnPromotersConsv.starts               = aTFBSOnPromotersConsv.pMotif, aTFBSOnPromotersConsv.pPromoters, aTFBSOnPromotersConsv.pConsvPerm
aTFBSOnPromotersConsv.ends                 = aTFBSOnPromotersConsv.pConsv
aTFBSOnPromotersConsv.pConsv.depends       = aTFBSOnPromotersConsv.pMotifScan, aTFBSOnPromotersConsv.pConsvPerm
aTFBSOnPromotersConsv.pMotifScan.depends   = aTFBSOnPromotersConsv.pMotif, aTFBSOnPromotersConsv.pBedGetfasta
aTFBSOnPromotersConsv.pBedGetfasta.depends = aTFBSOnPromotersConsv.pPromoters
# input
aTFBSOnPromotersConsv.pConsv.input = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTFBSOnPromotersConsv.pBedGetfasta.args.params.name = True
aTFBSOnPromotersConsv.pMotifScan.args.gnorm         = True
aTFBSOnPromotersConsv.pConsv.args.pval              = True


"""
@name:
	aTFBSOnRegionsConsv
@description:
	Scan motifs on a given regions with conservation score.
@depends:
	pBedGetfasta[*]  \
						pMotifScan
		   pMotif[*] /             \
		                            pConsv[!]
				             	   /
					  pConsvPerm[*]
"""
aTFBSOnRegionsConsv = Aggr(
	pFile2Proc.copy(newid = 'pMotif'),
	pBedGetfasta,
	pConsvPerm,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTFBSOnRegionsConsv.pMotif.runner     = 'local'
aTFBSOnRegionsConsv.pConsvPerm.input  = [0]
# delegate
aTFBSOnRegionsConsv.delegate('args.mids', 'pMotifScan')
aTFBSOnRegionsConsv.delegate('args.gnorm', 'pMotifScan')
aTFBSOnRegionsConsv.delegate('args.ref', 'pBedGetfasta')
aTFBSOnRegionsConsv.delegate('args.pval', 'pMotifScan')
aTFBSOnRegionsConsv.delegate('args.threspval', 'pConsv')
aTFBSOnRegionsConsv.delegate('args.len', 'pConsvPerm')
aTFBSOnRegionsConsv.delegate('args.nperm', 'pConsvPerm')
aTFBSOnRegionsConsv.delegate('args.chrsizes', 'pConsvPerm')
# depends
aTFBSOnRegionsConsv.starts             = aTFBSOnRegionsConsv.pMotif, aTFBSOnRegionsConsv.pBedGetfasta, aTFBSOnRegionsConsv.pConsvPerm
aTFBSOnRegionsConsv.ends               = aTFBSOnRegionsConsv.pConsv
aTFBSOnRegionsConsv.pConsv.depends     = aTFBSOnRegionsConsv.pMotifScan, aTFBSOnRegionsConsv.pConsvPerm
aTFBSOnRegionsConsv.pMotifScan.depends = aTFBSOnRegionsConsv.pMotif, aTFBSOnRegionsConsv.pBedGetfasta
# input
aTFBSOnRegionsConsv.pConsv.input = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTFBSOnRegionsConsv.pBedGetfasta.args.params.name = True
aTFBSOnRegionsConsv.pMotifScan.args.gnorm         = True
aTFBSOnRegionsConsv.pConsv.args.pval              = True



