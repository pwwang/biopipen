from pyppl import Aggr
from bioprocs.common import pStr2File, pFile2Proc
from bioprocs.seq import pPromoters
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



