from collections import defaultdict
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
tmbtype = {{args.type | quote}}

"""
Red --> Nonsense_Mutation, frameshift_variant, stop_gained, splice_acceptor_variant, splice_acceptor_variant&intron_variant, splice_donor_variant, splice_donor_variant&intron_variant, Splice_Site, Frame_Shift_Del, Frame_Shift_Ins

Blue --> splice_region_variant, splice_region_variant&intron_variant, missense, non_coding_exon_variant, missense_variant, Missense_Mutation, exon_variant, RNA, Indel, start_lost, start_gained, De_novo_Start_OutOfFrame, Translation_Start_Site, De_novo_Start_InFrame, stop_lost, Nonstop_Mutation, initiator_codon_variant, 5_prime_UTR_premature_start_codon_gain_variant, disruptive_inframe_deletion, inframe_deletion, inframe_insertion, In_Frame_Del, In_Frame_Ins

Green --> synonymous_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, 5'Flank, 3'Flank, 3'UTR, 5'UTR, Silent, stop_retained_variant

Orange --> others, SV, upstreamgenevariant, downstream_gene_variant, intron_variant, intergenic_region
"""
# oncotator MAFs
INTRON = "Intron"
FIVE_PRIME_UTR = "5'UTR"
THREE_PRIME_UTR = "3'UTR"
IGR = "IGR"
FIVE_PRIME_PRIME_FLANK = "5'Flank"
THREE_PRIME_PRIME_FLANK = "3'Flank"
MISSENSE = "Missense_Mutation"
NONSENSE = "Nonsense_Mutation"
NONSTOP = "Nonstop_Mutation"
SILENT = "Silent"
SPLICE_SITE = "Splice_Site"
IN_FRAME_DEL = "In_Frame_Del"
IN_FRAME_INS = "In_Frame_Ins"
FRAME_SHIFT_INS = "Frame_Shift_Ins"
FRAME_SHIFT_DEL = "Frame_Shift_Del"
START_CODON_SNP = "Start_Codon_SNP"
START_CODON_INS = "Start_Codon_Ins"
START_CODON_DEL = "Start_Codon_Del"
STOP_CODON_INS = "Stop_Codon_Ins"
STOP_CODON_DEL = "Stop_Codon_Del"
# Note: A STOP_CODON_SNP is a nonstop mutation (or Silent)
DE_NOVO_START_IN_FRAME = "De_novo_Start_InFrame"
DE_NOVO_START_OUT_FRAME = "De_novo_Start_OutOfFrame"
RNA = "RNA"
LINCRNA = "lincRNA"
# vcf2maf
TRANSLATION_START_SITE = 'Translation_Start_Site'
SPLICE_REGION = "Splice_Region"
TARGETED_REGION = "Targeted_Region"

def nonsyn(maf, outfile):
	# should only include Red and Blue
	excludes = (INTRON, FIVE_PRIME_UTR, THREE_PRIME_UTR, IGR, FIVE_PRIME_PRIME_FLANK, THREE_PRIME_PRIME_FLANK, SILENT, )
	reader = TsvReader(infile)
	# Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS      dbSNP_Val_Status        Tumor_Sample_Barcode ...
	tmb = defaultdict(lambda: 0)
	for r in reader:
		if r.Variant_Classification in excludes:
			continue
		tmb[r.Tumor_Sample_Barcode] += 1
	writer = TsvWriter(outfile)
	writer.cnames = ['Sample', 'TMB']
	writer.writeHead()
	for key, val in tmb.items():
		writer.write([key, val])

globals()[tmbtype](infile, outfile)