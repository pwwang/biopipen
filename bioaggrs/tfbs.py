from pyppl import Aggr
from bioprocs.common import pStr2File, pFile2Proc, pSort, pSimRead
from bioprocs.seq import pPromoters, pConsv, pConsvPerm
from bioprocs.bed import pBedGetfasta
from bioprocs.tfbs import pMotifScan
from bioprocs import params

"""
@name:
	aTFBSOnPromotersByTF
@description:
	Scan motifs on genes' promoter regions by giving TF names.
@depends:
	          pPromoters[*] \
	                        |
	                   pBedGetfasta \
	                                  pMotifScan[!]
	        pSortTFs[*] -- pSimRead /
	                          |
	pTFList[*] -- pSortTFList /
@input:
	- TF list file, one per line
	- gene list file, one per line
	- TF (1st col) list file with motif ids (2nd col). Default: params.tflist.value
"""
aTFBSOnPromotersByTF = Aggr(
	pSort.copy(newid = 'pSortTFs'),
	pPromoters,
	pFile2Proc.copy(newid = 'pTFList'),
	pSort.copy(newid = 'pSortTFList'),
	pSimRead,
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTFBSOnPromotersByTF.pTFList.runner    = 'local'
aTFBSOnPromotersByTF.pPromoters.runner = 'local'
# delegate
aTFBSOnPromotersByTF.delegate('args.up'      , 'pPromoters')
aTFBSOnPromotersByTF.delegate('args.down'    , 'pPromoters')
aTFBSOnPromotersByTF.delegate('args.genome'  , 'pPromoters')
aTFBSOnPromotersByTF.delegate('args.ref'     , 'pBedGetfasta')
aTFBSOnPromotersByTF.delegate('args.pval'    , 'pMotifScan')
aTFBSOnPromotersByTF.delegate('args.tfmotifs', 'pMotifScan')
# depends
aTFBSOnPromotersByTF.starts               = aTFBSOnPromotersByTF.pSortTFs, aTFBSOnPromotersByTF.pPromoters, aTFBSOnPromotersByTF.pTFList 
aTFBSOnPromotersByTF.ends                 = aTFBSOnPromotersByTF.pMotifScan
aTFBSOnPromotersByTF.pMotifScan.depends   = aTFBSOnPromotersByTF.pSimRead, aTFBSOnPromotersByTF.pBedGetfasta
aTFBSOnPromotersByTF.pBedGetfasta.depends = aTFBSOnPromotersByTF.pPromoters
aTFBSOnPromotersByTF.pSimRead.depends     = aTFBSOnPromotersByTF.pSortTFs, aTFBSOnPromotersByTF.pSortTFList
aTFBSOnPromotersByTF.pSortTFList.depends  = aTFBSOnPromotersByTF.pTFList
# input
aTFBSOnPromotersByTF.pTFList.input    = [params.tflist.value]
aTFBSOnPromotersByTF.pSimRead.input   = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)).flatten() for l in [max(ch1.length(), ch2.length())]]
aTFBSOnPromotersByTF.pMotifScan.input = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)) for l in [max(ch1.length(), ch2.length())]][0]
# args
aTFBSOnPromotersByTF.pBedGetfasta.args.params.name = True
aTFBSOnPromotersByTF.pSortTFList.args.params.k     = 2
aTFBSOnPromotersByTF.pSimRead.args.match           = 'lambda line1, line2: -1 if line1[0] == line2[1] else 0 if line1[0] < line2[1] else 1'
aTFBSOnPromotersByTF.pSimRead.args.do              = 'lambda line1, line2: fout.write("\\t".join(line2) + "\\n")'

"""
@name:
	aTFBSOnPromoters
@description:
	Scan motifs on genes' promoter regions.
@depends:
	pPromoters[*] \
					pBedGetfasta \
								  pMotifScan[!]
					  pTFList[*] /
@input:
	- TF (1st col) list file with motif ids (2nd col). Default: params.tflist.value
	- gene list file, one per line
"""
aTFBSOnPromoters = Aggr(
	pFile2Proc.copy(newid = 'pTFList'),
	pPromoters,
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTFBSOnPromoters.pTFList.runner    = 'local'
aTFBSOnPromoters.pPromoters.runner = 'local'
# delegate
aTFBSOnPromoters.delegate('args.up', 'pPromoters')
aTFBSOnPromoters.delegate('args.down', 'pPromoters')
aTFBSOnPromoters.delegate('args.genome', 'pPromoters')
aTFBSOnPromoters.delegate('args.ref', 'pBedGetfasta')
aTFBSOnPromoters.delegate('args.pval', 'pMotifScan')
aTFBSOnPromoters.delegate('args.tfmotifs', 'pMotifScan')
# depends
aTFBSOnPromoters.starts               = aTFBSOnPromoters.pTFList,     aTFBSOnPromoters.pPromoters 
aTFBSOnPromoters.ends                 = aTFBSOnPromoters.pMotifScan
aTFBSOnPromoters.pMotifScan.depends   = aTFBSOnPromoters.pTFList,     aTFBSOnPromoters.pBedGetfasta
aTFBSOnPromoters.pBedGetfasta.depends = aTFBSOnPromoters.pPromoters
# args
aTFBSOnPromoters.pBedGetfasta.args.params.name = True


"""
@name:
	aTFBSOnRegionsByTF
@description:
	Scan motifs on a given regions.
@depends:
	                      pBedGetfasta[*]  \
						                     pMotifScan[!]
		   	        pSortTFs[*] -- pSimRead /
	                          |
	pTFList[*] -- pSortTFList /
@input:
	- TF list file, one per line
	- region file in bed
	- TF (1st col) list file with motif ids (2nd col). Default: params.tflist.value
"""
aTFBSOnRegionsByTF = Aggr(
	pSort.copy(newid = 'pSortTFs'),
	pBedGetfasta,
	pFile2Proc.copy(newid = 'pTFList'),
	pSort.copy(newid = 'pSortTFList'),
	pSimRead,
	pMotifScan,
	depends = False
)
# defaults
aTFBSOnRegionsByTF.pTFList.runner     = 'local'
# delegate
aTFBSOnRegionsByTF.delegate('args.ref', 'pBedGetfasta')
aTFBSOnRegionsByTF.delegate('args.pval', 'pMotifScan')
# depends
aTFBSOnRegionsByTF.starts              = aTFBSOnRegionsByTF.pSortTFs, aTFBSOnRegionsByTF.pBedGetfasta, aTFBSOnRegionsByTF.pTFList
aTFBSOnRegionsByTF.ends                = aTFBSOnRegionsByTF.pMotifScan
aTFBSOnRegionsByTF.pMotifScan.depends  = aTFBSOnRegionsByTF.pSimRead, aTFBSOnRegionsByTF.pBedGetfasta
aTFBSOnRegionsByTF.pSimRead.depends    = aTFBSOnRegionsByTF.pSortTFs, aTFBSOnRegionsByTF.pSortTFList
aTFBSOnRegionsByTF.pSortTFList.depends = aTFBSOnRegionsByTF.pTFList
# input
aTFBSOnRegionsByTF.pTFList.input    = [params.tflist.value]
aTFBSOnRegionsByTF.pSimRead.input   = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)).flatten() for l in [max(ch1.length(), ch2.length())]]
aTFBSOnRegionsByTF.pMotifScan.input = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)) for l in [max(ch1.length(), ch2.length())]][0]
# args
aTFBSOnRegionsByTF.pBedGetfasta.args.params.name = True
aTFBSOnRegionsByTF.pSortTFList.args.params.k     = 2
aTFBSOnRegionsByTF.pSimRead.args.match           = 'lambda line1, line2: -1 if line1[0] == line2[1] else 0 if line1[0] < line2[1] else 1'
aTFBSOnRegionsByTF.pSimRead.args.do              = 'lambda line1, line2: fout.write("\\t".join(line2) + "\\n")'


"""
@name:
	aTFBSOnRegions
@description:
	Scan motifs on a given regions.
@depends:
	pBedGetfasta[*]  \
						pMotifScan[!]
		   pTFList[*] /
@input:
	- TF (1st col) list file with motif ids (2nd col). Default: params.tflist.value
	- region file in bed
"""
aTFBSOnRegions = Aggr(
	pFile2Proc.copy(newid = 'pTFList'),
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTFBSOnRegions.pTFList.runner     = 'local'
# delegate
aTFBSOnRegions.delegate('args.ref', 'pBedGetfasta')
aTFBSOnRegions.delegate('args.pval', 'pMotifScan')
# depends
aTFBSOnRegions.starts             = aTFBSOnRegions.pTFList, aTFBSOnRegions.pBedGetfasta
aTFBSOnRegions.ends               = aTFBSOnRegions.pMotifScan
aTFBSOnRegions.pMotifScan.depends = aTFBSOnRegions.pTFList, aTFBSOnRegions.pBedGetfasta
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
					   pTFList[*] /	         \
					                          pConsv[!]
											 /
						    	pConsvPerm[*]
@input:
	- TF list file, one per line
	- gene list file, one per line
	- Seeds for pConsvPerm. Default: [0]
"""
aTFBSOnPromotersConsv = Aggr(
	pFile2Proc.copy(newid = 'pTFList'),
	pPromoters,
	pConsvPerm,
	pBedGetfasta,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTFBSOnPromotersConsv.pTFList.runner    = 'local'
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
aTFBSOnPromotersConsv.starts               = aTFBSOnPromotersConsv.pTFList, aTFBSOnPromotersConsv.pPromoters, aTFBSOnPromotersConsv.pConsvPerm
aTFBSOnPromotersConsv.ends                 = aTFBSOnPromotersConsv.pConsv
aTFBSOnPromotersConsv.pConsv.depends       = aTFBSOnPromotersConsv.pMotifScan, aTFBSOnPromotersConsv.pConsvPerm
aTFBSOnPromotersConsv.pMotifScan.depends   = aTFBSOnPromotersConsv.pTFList, aTFBSOnPromotersConsv.pBedGetfasta
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
		   pTFList[*] /             \
		                            pConsv[!]
				             	   /
					  pConsvPerm[*]
@input:
	- TF list file, one per line
	- region list in bed.
	- Seeds for pConsvPerm. Default: [0]
"""
aTFBSOnRegionsConsv = Aggr(
	pFile2Proc.copy(newid = 'pTFList'),
	pBedGetfasta,
	pConsvPerm,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTFBSOnRegionsConsv.pTFList.runner     = 'local'
aTFBSOnRegionsConsv.pConsvPerm.input  = [0]
# delegate
aTFBSOnRegionsConsv.delegate('args.ref', 'pBedGetfasta')
aTFBSOnRegionsConsv.delegate('args.pval', 'pMotifScan')
aTFBSOnRegionsConsv.delegate('args.threspval', 'pConsv')
aTFBSOnRegionsConsv.delegate('args.len', 'pConsvPerm')
aTFBSOnRegionsConsv.delegate('args.nperm', 'pConsvPerm')
aTFBSOnRegionsConsv.delegate('args.chrsizes', 'pConsvPerm')
# depends
aTFBSOnRegionsConsv.starts             = aTFBSOnRegionsConsv.pTFList, aTFBSOnRegionsConsv.pBedGetfasta, aTFBSOnRegionsConsv.pConsvPerm
aTFBSOnRegionsConsv.ends               = aTFBSOnRegionsConsv.pConsv
aTFBSOnRegionsConsv.pConsv.depends     = aTFBSOnRegionsConsv.pMotifScan, aTFBSOnRegionsConsv.pConsvPerm
aTFBSOnRegionsConsv.pMotifScan.depends = aTFBSOnRegionsConsv.pTFList, aTFBSOnRegionsConsv.pBedGetfasta
# input
aTFBSOnRegionsConsv.pConsv.input = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTFBSOnRegionsConsv.pBedGetfasta.args.params.name = True
aTFBSOnRegionsConsv.pConsv.args.pval              = True



