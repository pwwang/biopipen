from pyppl import Aggr
from bioprocs.common import pStr2File, pFile2Proc, pSort
from bioprocs.matrix import pSimRead
from bioprocs.seq import pPromoters, pConsv, pConsvPerm
from bioprocs.bed import pBedGetfasta
from bioprocs.tfbs import pMotifScan
from bioprocs import params

"""
@name:
	aTfbsTfP
@description:
	Scan motifs on genes' promoter regions by giving TF names.
@depends:
		pPromoters[*] \
					|
				pBedGetfasta \
								pMotifScan[!]
	     pTFs[*] -- pSimRead /
						|
		    pTFList[*] /
@input:
	- TF list file, one per line
	- gene list file, one per line
	- TF (1st col) list file with motif ids (2nd col). Default: params.tflist.value
"""
aTfbsTfP = Aggr(
	pSort.copy(newid = 'pTFs'),
	pPromoters,
	pSort.copy(newid = 'pTFList'),
	pSimRead,
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTfbsTfP.pPromoters.runner = 'local'
# delegate
aTfbsTfP.delegate('args.up'      , 'pPromoters')
aTfbsTfP.delegate('args.down'    , 'pPromoters')
aTfbsTfP.delegate('args.genome'  , 'pPromoters')
aTfbsTfP.delegate('args.ref'     , 'pBedGetfasta')
aTfbsTfP.delegate('args.pval'    , 'pMotifScan')
aTfbsTfP.delegate('args.tfmotifs', 'pMotifScan')
# depends
aTfbsTfP.starts               = aTfbsTfP.pTFs, aTfbsTfP.pPromoters, aTfbsTfP.pTFList 
aTfbsTfP.ends                 = aTfbsTfP.pMotifScan
aTfbsTfP.pMotifScan.depends   = aTfbsTfP.pSimRead, aTfbsTfP.pBedGetfasta
aTfbsTfP.pBedGetfasta.depends = aTfbsTfP.pPromoters
aTfbsTfP.pSimRead.depends     = aTfbsTfP.pTFList, aTfbsTfP.pTFs
# input
aTfbsTfP.pTFList.input    = [params.tflist.value]
# input size of either pTFs or pTFList should be 1
aTfbsTfP.pSimRead.input   = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)).flatten() for l in [max(ch1.length(), ch2.length())]]
aTfbsTfP.pMotifScan.input = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)) for l in [max(ch1.length(), ch2.length())]][0]
# args
aTfbsTfP.pBedGetfasta.args.params.name = True
aTfbsTfP.pTFList.args.params.k         = 2
aTfbsTfP.pSimRead.args.match           = 'lambda line1, line2: -1 if line1[1] == line2[0] else 0 if line1[1] < line2[0] else 1'
aTfbsTfP.pSimRead.args.do              = 'lambda line1, line2: fout.write("\\t".join(line1) + "\\n")'

"""
@name:
	aTfbsP
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
aTfbsP = Aggr(
	pFile2Proc.copy(newid = 'pTFList'),
	pPromoters,
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTfbsP.pTFList.runner    = 'local'
aTfbsP.pPromoters.runner = 'local'
# delegate
aTfbsP.delegate('args.up', 'pPromoters')
aTfbsP.delegate('args.down', 'pPromoters')
aTfbsP.delegate('args.genome', 'pPromoters')
aTfbsP.delegate('args.ref', 'pBedGetfasta')
aTfbsP.delegate('args.pval', 'pMotifScan')
aTfbsP.delegate('args.tfmotifs', 'pMotifScan')
# depends
aTfbsP.starts               = aTfbsP.pTFList,     aTfbsP.pPromoters 
aTfbsP.ends                 = aTfbsP.pMotifScan
aTfbsP.pMotifScan.depends   = aTfbsP.pTFList,     aTfbsP.pBedGetfasta
aTfbsP.pBedGetfasta.depends = aTfbsP.pPromoters
# args
aTfbsP.pBedGetfasta.args.params.name = True


"""
@name:
	aTfbsTfR
@description:
	Scan motifs on a given regions.
@depends:
	                      pBedGetfasta[*]  \
						                     pMotifScan[!]
		   	           pTFs[*] -- pSimRead /
	                          |
	          pSortTFList[*] /
@input:
	- TF list file, one per line
	- region file in bed
	- TF (1st col) list file with motif ids (2nd col). Default: params.tflist.value
"""
aTfbsTfR = Aggr(
	pSort.copy(newid = 'pTFs'),
	pBedGetfasta,
	pSort.copy(newid = 'pTFList'),
	pSimRead,
	pMotifScan,
	depends = False
)
# defaults
aTfbsTfR.pTFList.runner     = 'local'
# delegate
aTfbsTfR.delegate('args.ref', 'pBedGetfasta')
aTfbsTfR.delegate('args.pval', 'pMotifScan')
aTfbsTfR.delegate('args.tfmotifs', 'pMotifScan')
# depends
aTfbsTfR.starts              = aTfbsTfR.pTFs, aTfbsTfR.pBedGetfasta, aTfbsTfR.pTFList
aTfbsTfR.ends                = aTfbsTfR.pMotifScan
aTfbsTfR.pMotifScan.depends  = aTfbsTfR.pSimRead, aTfbsTfR.pBedGetfasta
aTfbsTfR.pSimRead.depends    = aTfbsTfR.pTFList, aTfbsTfR.pTFs
# input
aTfbsTfR.pTFList.input    = [params.tflist.value]
aTfbsTfR.pSimRead.input   = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)).flatten() for l in [max(ch1.length(), ch2.length())]]
aTfbsTfR.pMotifScan.input = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)) for l in [max(ch1.length(), ch2.length())]][0]
# args
aTfbsTfR.pBedGetfasta.args.params.name = True
aTfbsTfR.pTFList.args.params.k         = 2
aTfbsTfR.pSimRead.args.match           = 'lambda line1, line2: -1 if line1[1] == line2[0] else 0 if line1[1] < line2[0] else 1'
aTfbsTfR.pSimRead.args.do              = 'lambda line1, line2: fout.write("\\t".join(line1) + "\\n")'


"""
@name:
	aTfbsR
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
aTfbsR = Aggr(
	pFile2Proc.copy(newid = 'pTFList'),
	pBedGetfasta,
	pMotifScan,
	depends = False
)
# defaults
aTfbsR.pTFList.runner     = 'local'
# delegate
aTfbsR.delegate('args.ref', 'pBedGetfasta')
aTfbsR.delegate('args.pval', 'pMotifScan')
aTfbsR.delegate('args.tfmotifs', 'pMotifScan')
# depends
aTfbsR.starts             = aTfbsR.pTFList, aTfbsR.pBedGetfasta
aTfbsR.ends               = aTfbsR.pMotifScan
aTfbsR.pMotifScan.depends = aTfbsR.pTFList, aTfbsR.pBedGetfasta
# args
aTfbsR.pBedGetfasta.args.params.name = True

"""
@name:
	aTfbsPC
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
aTfbsPC = Aggr(
	pFile2Proc.copy(newid = 'pTFList'),
	pPromoters,
	pConsvPerm,
	pBedGetfasta,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTfbsPC.pTFList.runner    = 'local'
aTfbsPC.pPromoters.runner = 'local'
aTfbsPC.pConsvPerm.input  = [0]
# delegate
aTfbsPC.delegate('args.up', 'pPromoters')
aTfbsPC.delegate('args.down', 'pPromoters')
aTfbsPC.delegate('args.genome', 'pPromoters')
aTfbsPC.delegate('args.ref', 'pBedGetfasta')
aTfbsPC.delegate('args.pval', 'pMotifScan')
aTfbsPC.delegate('args.tfmotifs', 'pMotifScan')
aTfbsPC.delegate('args.cpval', 'pConsv', 'args.pval')
aTfbsPC.delegate('args.len', 'pConsvPerm')
aTfbsPC.delegate('args.nperm', 'pConsvPerm')
aTfbsPC.delegate('args.chrsizes', 'pConsvPerm')
# depends
aTfbsPC.starts               = aTfbsPC.pTFList, aTfbsPC.pPromoters, aTfbsPC.pConsvPerm
aTfbsPC.ends                 = aTfbsPC.pConsv
aTfbsPC.pConsv.depends       = aTfbsPC.pMotifScan, aTfbsPC.pConsvPerm
aTfbsPC.pMotifScan.depends   = aTfbsPC.pTFList, aTfbsPC.pBedGetfasta
aTfbsPC.pBedGetfasta.depends = aTfbsPC.pPromoters
# input
aTfbsPC.pConsv.input = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTfbsPC.pBedGetfasta.args.params.name = True
aTfbsPC.pConsv.args.pval              = 0.05


"""
@name:
	aTfbsRC
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
aTfbsRC = Aggr(
	pFile2Proc.copy(newid = 'pTFList'),
	pBedGetfasta,
	pConsvPerm,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTfbsRC.pTFList.runner     = 'local'
aTfbsRC.pConsvPerm.input  = [0]
# delegate
aTfbsRC.delegate('args.ref', 'pBedGetfasta')
aTfbsRC.delegate('args.pval', 'pMotifScan')
aTfbsRC.delegate('args.tfmotifs', 'pMotifScan')
aTfbsRC.delegate('args.cpval', 'pConsv', 'args.pval')
aTfbsRC.delegate('args.len', 'pConsvPerm')
aTfbsRC.delegate('args.nperm', 'pConsvPerm')
aTfbsRC.delegate('args.chrsizes', 'pConsvPerm')
# depends
aTfbsRC.starts             = aTfbsRC.pTFList, aTfbsRC.pBedGetfasta, aTfbsRC.pConsvPerm
aTfbsRC.ends               = aTfbsRC.pConsv
aTfbsRC.pConsv.depends     = aTfbsRC.pMotifScan, aTfbsRC.pConsvPerm
aTfbsRC.pMotifScan.depends = aTfbsRC.pTFList, aTfbsRC.pBedGetfasta
# input
aTfbsRC.pConsv.input = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTfbsRC.pBedGetfasta.args.params.name = True
aTfbsRC.pConsv.args.pval              = 0.05

"""
@name:
	aTfbsTfPC
@description:
	Scan motifs on genes' promoter regions by giving TF names with conservation.
@depends:
		pPromoters[*] \
					|
				pBedGetfasta \
								pMotifScan
	    pTFs[*] -- pSimRead /               \
						|                    \
		    pTFList[*] /                       pConsv[!]
			                   pConsvPerm[*] /
@input:
	- TF list file, one per line
	- gene list file, one per line
	- TF (1st col) list file with motif ids (2nd col). Default: params.tflist.value
"""
aTfbsTfPC = Aggr(
	pSort.copy(newid = 'pTFs'),
	pPromoters,
	pSort.copy(newid = 'pTFList'),
	pSimRead,
	pBedGetfasta,
	pConsvPerm,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTfbsTfPC.pPromoters.runner = 'local'
aTfbsTfPC.pConsvPerm.input  = [0]
# delegate
aTfbsTfPC.delegate('args.up'      , 'pPromoters')
aTfbsTfPC.delegate('args.down'    , 'pPromoters')
aTfbsTfPC.delegate('args.genome'  , 'pPromoters')
aTfbsTfPC.delegate('args.ref'     , 'pBedGetfasta')
aTfbsTfPC.delegate('args.pval'    , 'pMotifScan')
aTfbsTfPC.delegate('args.tfmotifs', 'pMotifScan')
aTfbsTfPC.delegate('args.cpval'   , 'pConsv', 'args.pval')
aTfbsTfPC.delegate('args.len'     , 'pConsvPerm')
aTfbsTfPC.delegate('args.nperm'   , 'pConsvPerm')
aTfbsTfPC.delegate('args.chrsizes', 'pConsvPerm')
# depends
aTfbsTfPC.starts               = aTfbsTfPC.pTFs, aTfbsTfPC.pPromoters, aTfbsTfPC.pTFList, aTfbsTfPC.pConsvPerm 
aTfbsTfPC.ends                 = aTfbsTfPC.pConsv
aTfbsTfPC.pConsv.depends       = aTfbsTfPC.pMotifScan, aTfbsTfPC.pConsvPerm
aTfbsTfPC.pMotifScan.depends   = aTfbsTfPC.pSimRead, aTfbsTfPC.pBedGetfasta
aTfbsTfPC.pBedGetfasta.depends = aTfbsTfPC.pPromoters
aTfbsTfPC.pSimRead.depends     = aTfbsTfPC.pTFList, aTfbsTfPC.pTFs
# input
aTfbsTfPC.pTFList.input    = [params.tflist.value]
# input size of either pTFs or pTFList should be 1
aTfbsTfPC.pSimRead.input   = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)).flatten() for l in [max(ch1.length(), ch2.length())]]
aTfbsTfPC.pMotifScan.input = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)) for l in [max(ch1.length(), ch2.length())]][0]
aTfbsTfPC.pConsv.input     = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTfbsTfPC.pBedGetfasta.args.params.name = True
aTfbsTfPC.pTFList.args.params.k         = 2
aTfbsTfPC.pSimRead.args.match           = 'lambda line1, line2: -1 if line1[1] == line2[0] else 0 if line1[1] < line2[0] else 1'
aTfbsTfPC.pSimRead.args.do              = 'lambda line1, line2: fout.write("\\t".join(line1) + "\\n")'
aTfbsTfPC.pConsv.args.pval              = 0.05

"""
@name:
	aTfbsTfRC
@description:
	Scan motifs on a given regions with conservation.
@depends:
	                      pBedGetfasta[*]  \
						                     pMotifScan
		   	            pTFs[*] -- pSimRead /         \
	                          |                        \
	          pSortTFList[*] /                           pConsv[!]
			                            pConsvPerm[*]  /
@input:
	- TF list file, one per line
	- region file in bed
	- TF (1st col) list file with motif ids (2nd col). Default: params.tflist.value
"""
aTfbsTfRC = Aggr(
	pSort.copy(newid = 'pTFs'),
	pBedGetfasta,
	pConsvPerm,
	pSort.copy(newid = 'pTFList'),
	pSimRead,
	pMotifScan,
	pConsv,
	depends = False
)
# defaults
aTfbsTfRC.pTFList.runner    = 'local'
aTfbsTfRC.pConsvPerm.input  = [0]
# delegate
aTfbsTfRC.delegate('args.ref', 'pBedGetfasta')
aTfbsTfRC.delegate('args.pval', 'pMotifScan')
aTfbsTfRC.delegate('args.tfmotifs', 'pMotifScan')
aTfbsTfRC.delegate('args.cpval'   , 'pConsv', 'args.pval')
aTfbsTfRC.delegate('args.len'     , 'pConsvPerm')
aTfbsTfRC.delegate('args.nperm'   , 'pConsvPerm')
aTfbsTfRC.delegate('args.chrsizes', 'pConsvPerm')
# depends
aTfbsTfRC.starts             = aTfbsTfRC.pTFs,   aTfbsTfRC.pBedGetfasta, aTfbsTfRC.pTFList, aTfbsTfRC.pConsvPerm
aTfbsTfRC.ends               = aTfbsTfRC.pConsv
aTfbsTfRC.pConsv.depends     = aTfbsTfRC.pMotifScan, aTfbsTfRC.pConsvPerm
aTfbsTfRC.pMotifScan.depends = aTfbsTfRC.pSimRead,   aTfbsTfRC.pBedGetfasta
aTfbsTfRC.pSimRead.depends   = aTfbsTfRC.pTFList,    aTfbsTfRC.pTFs
# input
aTfbsTfRC.pTFList.input    = [params.tflist.value]
aTfbsTfRC.pSimRead.input   = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)).flatten() for l in [max(ch1.length(), ch2.length())]]
aTfbsTfRC.pMotifScan.input = lambda ch1, ch2: [ch1.repRow(l).cbind(ch2.repRow(l)) for l in [max(ch1.length(), ch2.length())]][0]
aTfbsTfRC.pConsv.input     = lambda ch1, ch2: ch1.outfile.cbind(ch2)
# args
aTfbsTfRC.pBedGetfasta.args.params.name = True
aTfbsTfRC.pTFList.args.params.k         = 2
aTfbsTfRC.pSimRead.args.match           = 'lambda line1, line2: -1 if line1[1] == line2[0] else 0 if line1[1] < line2[0] else 1'
aTfbsTfRC.pSimRead.args.do              = 'lambda line1, line2: fout.write("\\t".join(line1) + "\\n")'
aTfbsTfRC.pConsv.args.pval              = 0.05

