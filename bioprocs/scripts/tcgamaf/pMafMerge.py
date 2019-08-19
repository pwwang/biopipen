from os import path
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

maffiles = {{i.infiles | repr}}
excols   = {{args.excols | quote}}
outfile  = {{o.outfile | quote}}

# first 34 columns are TCGA standard columns
headers = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "Tumor_Validation_Allele1", "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2", "Verification_Status", "Validation_Status", "Mutation_Status", "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score", "BAM_File", "Sequencer", "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID"]

# get all headers
cols   = []
for maffile in maffiles:
	reader = TsvReader(maffile)
	headers.extend(cname for cname in reader.cnames if cname not in headers)
	reader.close()

if excols == 'discard':
	headers = headers[:34]

writer = TsvWriter(outfile)
writer.cnames = headers
writer.writeHead()
for maffile in maffiles:
	reader = TsvReader(maffile)
	for r in reader:
		for head in headers:
			if head not in r:
				r[head] = '__UNKNOWN__'
		writer.write(r)
	reader.close()
writer.close()