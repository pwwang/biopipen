from os import path
from pyppl.utils import always_list
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile    = {{i.infile | quote}}
classfile = {{i.classfile | quote}}
outfile   = {{o.outfile | quote}}
classes   = {{args.classes | repr}}
"""
Red --> Nonsense_Mutation, frameshift_variant, stop_gained, splice_acceptor_variant, splice_acceptor_variant&intron_variant, splice_donor_variant, splice_donor_variant&intron_variant, Splice_Site, Frame_Shift_Del, Frame_Shift_Ins

Blue --> splice_region_variant, splice_region_variant&intron_variant, missense, non_coding_exon_variant, missense_variant, Missense_Mutation, exon_variant, RNA, Indel, start_lost, start_gained, De_novo_Start_OutOfFrame, Translation_Start_Site, De_novo_Start_InFrame, stop_lost, Nonstop_Mutation, initiator_codon_variant, 5_prime_UTR_premature_start_codon_gain_variant, disruptive_inframe_deletion, inframe_deletion, inframe_insertion, In_Frame_Del, In_Frame_Ins

Green --> synonymous_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, 5'Flank, 3'Flank, 3'UTR, 5'UTR, Silent, stop_retained_variant

Orange --> others, SV, upstreamgenevariant, downstream_gene_variant, intron_variant, intergenic_region
"""

CLASSES = dict(
	# oncotator MAFs
	INTRON                  = "Intron", # Orange
	FIVE_PRIME_UTR          = "5'UTR", # Green
	THREE_PRIME_UTR         = "3'UTR", # Green
	IGR                     = "IGR", # Orange
	FIVE_PRIME_PRIME_FLANK  = "5'Flank", # Green
	THREE_PRIME_PRIME_FLANK = "3'Flank", # Green
	MISSENSE                = "Missense_Mutation", # Blue
	NONSENSE                = "Nonsense_Mutation", # Red
	NONSTOP                 = "Nonstop_Mutation", # Blue
	SILENT                  = "Silent", # Green
	SPLICE_SITE             = "Splice_Site", # Red
	IN_FRAME_DEL            = "In_Frame_Del", # Blue
	IN_FRAME_INS            = "In_Frame_Ins", # Blue
	FRAME_SHIFT_INS         = "Frame_Shift_Ins", # Red
	FRAME_SHIFT_DEL         = "Frame_Shift_Del", # Red
	START_CODON_SNP         = "Start_Codon_SNP", # Blue
	START_CODON_INS         = "Start_Codon_Ins", # Blue
	START_CODON_DEL         = "Start_Codon_Del", # Blue
	STOP_CODON_INS          = "Stop_Codon_Ins", # Blue
	STOP_CODON_DEL          = "Stop_Codon_Del", # Blue
	# Note: A STOP_CODON_SNP is a nonstop mutation (or Silent)
	DE_NOVO_START_IN_FRAME  = "De_novo_Start_InFrame", # Blue
	DE_NOVO_START_OUT_FRAME = "De_novo_Start_OutOfFrame", # Blue
	RNA                     = "RNA", # Blue
	LINCRNA                 = "lincRNA", # Blue
	# vcf2maf
	TRANSLATION_START_SITE = 'Translation_Start_Site', # Blue
	SPLICE_REGION          = "Splice_Region", # Red
	TARGETED_REGION        = "Targeted_Region", # Blue
)

if classfile and path.isfile(classfile):
	with open(classfile) as f:
		classes = f.read().strip().splitlines()
elif classfile:
	classes = classfile.split(',')
else:
	classes = always_list(classes)

if any(c.upper() == 'RED' for c in classes):
	classes = [c for c in classes if c.upper() != 'RED']
	classes.extend([
		CLASSES['NONSENSE'],
		CLASSES['SPLICE_SITE'],
		CLASSES['FRAME_SHIFT_INS'],
		CLASSES['FRAME_SHIFT_DEL'],
		CLASSES['SPLICE_REGION']
	])

if any(c.upper() == 'BLUE' for c in classes):
	classes = [c for c in classes if c.upper() != 'BLUE']
	classes.extend([
		CLASSES['MISSENSE'],
		CLASSES['NONSTOP'],
		CLASSES['IN_FRAME_DEL'],
		CLASSES['IN_FRAME_INS'],
		CLASSES['START_CODON_SNP'],
		CLASSES['START_CODON_INS'],
		CLASSES['START_CODON_DEL'],
		CLASSES['STOP_CODON_INS'],
		CLASSES['STOP_CODON_DEL'],
		CLASSES['DE_NOVO_START_IN_FRAME'],
		CLASSES['DE_NOVO_START_OUT_FRAME'],
		CLASSES['RNA'],
		CLASSES['LINCRNA'],
		CLASSES['TRANSLATION_START_SITE'],
		CLASSES['TARGETED_REGION'],
	])

if any(c.upper() == 'GREEN' for c in classes):
	classes = [c for c in classes if c.upper() != 'GREEN']
	classes.extend([
		CLASSES['FIVE_PRIME_UTR'],
		CLASSES['THREE_PRIME_UTR'],
		CLASSES['FIVE_PRIME_PRIME_FLANK'],
		CLASSES['THREE_PRIME_PRIME_FLANK'],
		CLASSES['SILENT']
	])

if any(c.upper() == 'ORANGE' for c in classes):
	classes = [c for c in classes if c.upper() != 'ORANGE']
	classes.extend([
		CLASSES['INTRON'],
		CLASSES['IGR'],
	])

reader = TsvReader(infile)
writer = TsvWriter(outfile)
writer.cnames = reader.cnames
writer.writeHead()

for r in reader:
	if r.Variant_Classification in classes:
		writer.write(r)

reader.close()
writer.close()
