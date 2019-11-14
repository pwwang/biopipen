from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell, logger
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
dtoxog  = {{args.dtoxog | quote}}
params  = {{args.params | repr}}
keep    = {{args.keep | repr}}
nthread = {{args.nthread | repr}}

shell.load_config(dtoxog = dtoxog)

# see if maf file has required fields
# required fields
# Chromosome -- Contig number without a prefix.  E.g. "3" or "X" (without quotes)
# Start_position -- Position of the SNV
# End_position -- Should be the same as Start_position for SNV
# Reference_Allele -- The allele found in the reference genome at this position.  Single character representing the base, e.g. "G".
# Tumor_Seq_Allele1 -- alternate allele.  Single character representing the base, e.g. "G".
# Tumor_Sample_Barcode -- name of the tumor (or "case") sample.  This is used to generate file names and plots.
# Matched_Norm_Sample_Barcode -- name of the normal (or "control") sample.  This is used to generate file names and plots.
# ref_context -- Small window into the reference at the SNV.  The center position should be the same as Reference_Allele.  The total string should be of odd length and have a minimum length of 3.
# 	For example: Reference Allele is G, Chromosome is 1, Start_position and End_position are 120906037:  ref_context is CTTTTTTCGCGCAAAAATGCC  (string size is 21, in this case)
# i_t_ALT_F1R2 -- the number of reads with pair orientation of F1R2 and with the alternate allele (Tumor_Seq_Allele1).
# i_t_ALT_F2R1 -- the number of reads with pair orientation of F2R1 and with the alternate allele (Tumor_Seq_Allele1).
# i_t_REF_F1R2 -- the number of reads with pair orientation of F1R2 and with the reference allele (Reference_Allele).
# i_t_REF_F2R1 -- the number of reads with pair orientation of F2R1 and with the reference allele (Reference_Allele).
# i_t_Foxog -- Foxog, as described in the methods.  Depending on the nature of the reference and alternate alleles, either i_t_ALT_F1R2/(i_t_ALT_F1R2 + i_t_ALT_F2R1)  or i_t_ALT_F2R1/(i_t_ALT_F1R2 + i_t_ALT_F2R1).
# 	C>anything:  numerator is i_t_ALT_F2R1
# 	A>anything:  numerator is i_t_ALT_F2R1
# 	G>anything:  numerator is i_t_ALT_F1R2
# 	T>anything:  numerator is i_t_ALT_F1R2
# Variant_Type -- "SNP" (without quotes)
# i_picard_oxoQ -- will be 0

required_fields = { # set
	'Chromosome',
	'Start_Position',
	'End_Position',
	'Reference_Allele',
	'Tumor_Seq_Allele1',
	'Tumor_Sample_Barcode',
	'Matched_Norm_Sample_Barcode',
	'ref_context',
	'i_t_ALT_F1R2',
	'i_t_ALT_F2R1',
	'i_t_REF_F1R2',
	'i_t_REF_F2R1',
	'i_t_Foxog',
	'Variant_Type',
	'i_picard_oxoQ',
}

# check
logger.info('Checking required fields ...')
header = None
with open(infile) as f:
	for line in f:
		if line.startswith('Hugo_Symbol'):
			header = line.rstrip('\n').split('\t')
			break

if not header:
	raise ValueError('No header found for input MAF file.')

fields_not_found = required_fields - set(header)
if fields_not_found and fields_not_found != {'i_picard_oxoQ'}:
	raise ValueError('Required fields not found: %s' % (fields_not_found))


logger.info('Cleaning up input MAF file ...')
# Add i_picard_oxoQ as 0 if column does not exist
# TODO: add real i_picard_oxoQ
# Set i_t_ALT_F1R2, i_t_ALT_F2R1, i_t_REF_F1R2, i_t_REF_F2R1 as 0 if they are not reported
# Exclude mutations other than [SNP, DNP, TNP]
infile_fixed = path.join(path.dirname(outfile), 'fixed.' + path.basename(infile))
reader = TsvReader(infile, comment = '#', cnames = True)
writer = TsvWriter(infile_fixed)
if fields_not_found:
	writer.cnames = reader.cnames + ['i_picard_oxoQ']
writer.writeHead(lambda cnames: [
	'Start_position' if cname == 'Start_Position' else
	'End_position' if cname == 'End_Position' else cname
	for cname in cnames])
for r in reader:
	if r.Variant_Type not in ('SNP', 'DNP', 'TNP'):
		continue
	r.i_picard_oxoQ = 0
	for key in ('i_t_ALT_F1R2', 'i_t_ALT_F2R1', 'i_t_REF_F1R2', 'i_t_REF_F2R1'):
		if not r[key].isdigit():
			r[key] = 0
	writer.write(r)
writer.close()

params.mafFilename       = infile_fixed
params.outputMAFFilename = path.basename(outfile) + '.all'
params.outputDir         = path.dirname(outfile)
params.nthread           = nthread

# avoid lost values
for key, value in params.items():
	if isinstance(value, bool):
		params[key] = int(value)

shell.mkdir(path.join(path.dirname(outfile), 'figures'))
logger.info('Running DToxoG ...')
shell.fg.dtoxog(**params)


logger.info('Fixing output file format issues ...')
if not keep:
	# remove variants with isArtifactMode = 1
	with open(outfile + '.all') as fin, open(outfile, 'w') as fout:
		for line in fin:
			if line.startswith('#'):
				fout.write(line)
			else:
				break
	reader = TsvReader(outfile + '.all')
	writer = TsvWriter(outfile, append = True)
	writer.cnames = reader.cnames
	writer.writeHead(lambda cnames: [
		'Start_Position' if cname == 'Start_position' else
		'End_Position' if cname == 'End_position' else cname
		for cname in cnames])
	for r in reader:
		if r.isArtifactMode == '1':
			continue
		writer.write(r)
	reader.close()
	writer.close()
else:
	# get Start_Position and End_Position back
	with open(outfile + '.all') as fin, open(outfile, 'w') as fout:
		for line in fin:
			if line.startswith('Hugo_Symbol'):
				fout.write(line.replace('Start_position', 'Start_Position').replace('End_position', 'End_Position'))
			else:
				fout.write(line)
