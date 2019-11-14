from os import path
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile    = {{i.infile | quote}}
classfile = {{i.classfile | quote}}
outfile   = {{o.outfile | quote}}
classes   = {{args.classes | repr}}

if classfile and path.isfile(classfile):
	with open(classfile) as f:
		classes = f.read().strip().splitlines()
elif classfile:
	classes = classfile.split(',')
elif not isinstance(classes, list):
	classes = classes.split(',')

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

standard_classes = [INTRON, FIVE_PRIME_UTR, THREE_PRIME_UTR, IGR, FIVE_PRIME_PRIME_FLANK, THREE_PRIME_PRIME_FLANK, MISSENSE, NONSENSE, NONSTOP, SILENT, SPLICE_SITE, IN_FRAME_DEL, IN_FRAME_INS, FRAME_SHIFT_INS, FRAME_SHIFT_DEL, START_CODON_SNP, START_CODON_INS, START_CODON_DEL, STOP_CODON_INS, STOP_CODON_DEL, DE_NOVO_START_IN_FRAME, DE_NOVO_START_OUT_FRAME, RNA, LINCRNA, TRANSLATION_START_SITE, SPLICE_REGION, TARGETED_REGION]

if any(klass not in standard_classes for klass in classes):
	raise TypeError('Not a standard variant class. Expect one of %s' % standard_classes)

reader = TsvReader(infile)
writer = TsvWriter(outfile)
writer.cnames = reader.cnames
writer.writeHead()

for r in reader:
	if r.Variant_Classification in classes:
		writer.write(r)

reader.close()
writer.close()
