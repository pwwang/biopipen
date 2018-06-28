from os import path

# first 34 columns are TCGA standard columns
headers = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "Tumor_Validation_Allele1", "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2", "Verification_Status", "Validation_Status", "Mutation_Status", "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score", "BAM_File", "Sequencer", "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID"]
maffiles = {{in.infiles | repr}}

# get all headers
mafver = ''
cols   = []
for maffile in maffiles:
	with open(maffile) as f:
		for line in f:
			if line.startswith('#version') and not mafver:
				mafver = line # have \n
			if line.startswith('Hugo_Symbol'):
				header = line.strip().split('\t')
				cols.append(header)
				excols = [col for col in header if col not in headers]
				headers += excols

{% if args.excols | lambda x: x == 'discard' %}
headers = headers[:34]
{% endif %}

with open({{out.outfile | quote}}, 'w') as fout:
	fout.write(mafver)
	fout.write('\t'.join(headers) + '\n')
	for i, maffile in enumerate(maffiles):
		header = cols[i]
		with open(maffile) as f:
			for line in f:
				if line.startswith('#') or line.startswith('Hugo_Symbol'): continue
				parts = line.strip('\r\n').split('\t')
				parts = [parts[header.index(col)] if col in header else '' for k,col in enumerate(headers)]
				fout.write('\t'.join(parts) + '\n')


