from gff import Gff # python-gff
from diot import Diot
from pathlib import Path
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile   = Path({{i.infile | quote}})
outfile  = Path({{o.outfile | quote}})
genome   = {{args.genome | quote}}
missing  = {{args.missing | quote}}
tumor    = {{args.tumor | quote}}
normal   = {{args.normal | quote}}
inopts   = {{args.inopts | repr}}
full     = {{args.full | bool}}
refall   = {{args.refall | quote}}
bedtools = {{args.bedtools | quote}}
shell.load_config(bedtools = bedtools)

FULL_COLUMNS = [
	'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome',
	'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type',
	'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status',
	'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1',
	'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2',
	'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status',
	'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source',
	'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID',
	'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID',
	'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count',
	'n_alt_count', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type',
	'One_Consequence', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position',
	'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE',
	'TRANSCRIPT_STRAND', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL',
	'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen',
	'EXON', 'INTRON', 'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EAS_MAF',
	'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME',
	'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS',
	'TSL', 'HGVS_OFFSET', 'PHENO', 'MINIMISED', 'ExAC_AF', 'ExAC_AF_Adj', 'ExAC_AF_AFR',
	'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS',
	'GENE_PHENO', 'FILTER', 'CONTEXT', 'src_vcf_id', 'tumor_bam_uuid', 'normal_bam_uuid',
	'case_id', 'GDC_FILTER', 'COSMIC', 'MC3_Overlap', 'GDC_Validation_Status',
	'GDC_Valid_Somatic', 'vcf_region', 'vcf_info', 'vcf_format', 'vcf_tumor_gt', 'vcf_normal_gt'
]

ESSENTIAL_COLUMNS = FULL_COLUMNS[4:6] + [FULL_COLUMNS[10], FULL_COLUMNS[12]]

reader = TsvReader(infile, **inopts)
for ecol in ESSENTIAL_COLUMNS:
	if ecol not in reader.cnames:
		raise ValueError('Essential column {!r} not found.'.format(ecol))

gene_mappings = {}
if 'Hugo_Symbol' not in reader.cnames and 'Gene' not in reader.cnames:

	refgene = outfile.parent.joinpath('.refgene.bed')
	with open(refgene, 'w') as fout:
		for gff in Gff(refall):
			if gff['type'] != 'gene':
				continue
			fout.write(f"{gff['seqid']}\t{gff['start']}\t{gff['end']}\t{gff['attributes']['gene_name']}\t{gff['attributes']['gene_id']}")

	bedfile = outfile.parent.joinpath(infile.with_suffix('.bed'))
	with open(bedfile, 'w') as fbed:
		for r in reader:
			r.End_Position = r.get('End_Position', int(r.Start_Position) + 1)
			fbed.write([r.Chromosome, r.Start_Position, r.End_Position])

	res = shell.bedtools.intersect(a = bedfile, b = refgene, wa = True, wb = True).iter()
	for line in res:
		parts = line.strip().split('\t')
		gene_mappings['\t'.join(parts[:3])] = tuple(parts[:-2])

	reader.rewind()

writer = TsvWriter(outfile)
writer.cnames = FULL_COLUMNS[:32] if not full else FULL_COLUMNS
writer.writeHead()
for r in reader:
	r.End_Position = r.get('End_Position', int(r.Start_Position) + 1)
	genes = gene_mappings.get('\t'.join([r.Chromosome, str(r.Start_Position), str(r.End_Position)]))
	r.Hugo_Symbol = r.get('Hugo_Symbol', genes[0] if genes else missing)
	r.Entrez_Gene_Id = r.get('Entrez_Gene_Id', 0)
	r.Center = r.get('Center', missing)
	r.NCBI_Build = r.get('NCBI_Build', genome)
	r.Strand = r.get('Strand', '+')
	r.Variant_Classification = r.get('Variant_Classification', missing)
	r.Variant_Type = r.get('Variant_Type', missing)
	r.Tumor_Seq_Allele1 = r.Reference_Allele
	r.dbSNP_RS = r.get('dbSNP_RS', missing)
	r.dbSNP_Val_Status = r.get('dbSNP_Val_Status', missing)
	r.Tumor_Sample_Barcode = r.get('Tumor_Sample_Barcode', tumor)
	r.Matched_Norm_Sample_Barcode = r.get('Matched_Norm_Sample_Barcode', normal)
	r.Match_Norm_Seq_Allele1 = r.get('Match_Norm_Seq_Allele1', r.Reference_Allele)
	r.Match_Norm_Seq_Allele2 = r.get('Match_Norm_Seq_Allele2', r.Reference_Allele)
	r.Tumor_Validation_Allele1 = r.get('Tumor_Validation_Allele1', missing)
	r.Tumor_Validation_Allele2 = r.get('Tumor_Validation_Allele2', missing)
	r.Match_Norm_Validation_Allele1 = r.get('Match_Norm_Validation_Allele1', missing)
	r.Match_Norm_Validation_Allele2 = r.get('Match_Norm_Validation_Allele2', missing)
	r.Verification_Status = r.get('Verification_Status', missing)
	r.Validation_Status = r.get('Validation_Status', missing)
	r.Mutation_Status = r.get('Mutation_Status', missing)
	r.Sequencing_Phase = r.get('Sequencing_Phase', missing)
	r.Sequence_Source = r.get('Sequence_Source', missing)
	r.Validation_Method = r.get('Validation_Method', missing)
	r.Score = r.get('Score', missing)
	r.BAM_File = r.get('BAM_File', missing)
	r.Sequencer = r.get('Sequencer', missing)
	r.Tumor_Sample_UUID = r.get('Tumor_Sample_UUID', missing)
	r.Matched_Norm_Sample_UUID = r.get('Matched_Norm_Sample_UUID', missing)
	r.HGVSc = r.get('HGVSc', missing)
	r.HGVSp = r.get('HGVSp', missing)
	r.HGVSp_Short = r.get('HGVSp_Short', missing)
	r.Transcript_ID = r.get('Transcript_ID', missing)
	r.Exon_Number = r.get('Exon_Number', missing)
	r.t_depth = r.get('t_depth', missing)
	r.t_ref_count = r.get('t_ref_count', missing)
	r.t_alt_count = r.get('t_alt_count', missing)
	r.n_depth = r.get('n_depth', missing)
	r.n_ref_count = r.get('n_ref_count', missing)
	r.n_alt_count = r.get('n_alt_count', missing)
	r.all_effects = r.get('all_effects', missing)
	r.Allele = r.get('Allele', missing)
	r.Gene = r.get('Gene', r.Hugo_Symbol)
	r.Feature = r.get('Feature', missing)
	r.Feature_type = r.get('Feature_type', missing)
	r.One_Consequence = r.get('One_Consequence', missing)
	r.Consequence = r.get('Consequence', missing)
	r.cDNA_position = r.get('cDNA_position', missing)
	r.CDS_position = r.get('CDS_position', missing)
	r.Protein_position = r.get('Protein_position', missing)
	r.Amino_acids = r.get('Amino_acids', missing)
	r.Codons = r.get('Codons', missing)
	r.Existing_variation = r.get('Existing_variation', missing)
	r.ALLELE_NUM = r.get('ALLELE_NUM', missing)
	r.DISTANCE = r.get('DISTANCE', missing)
	r.TRANSCRIPT_STRAND = r.get('TRANSCRIPT_STRAND', missing)
	r.SYMBOL = r.get('SYMBOL', missing)
	r.SYMBOL_SOURCE = r.get('SYMBOL_SOURCE', missing)
	r.HGNC_ID = r.get('HGNC_ID', missing)
	r.BIOTYPE = r.get('BIOTYPE', missing)
	r.CANONICAL = r.get('CANONICAL', missing)
	r.CCDS = r.get('CCDS', missing)
	r.ENSP = r.get('ENSP', missing)
	r.SWISSPROT = r.get('SWISSPROT', missing)
	r.TREMBL = r.get('TREMBL', missing)
	r.UNIPARC = r.get('UNIPARC', missing)
	r.RefSeq = r.get('RefSeq', missing)
	r.SIFT = r.get('SIFT', missing)
	r.PolyPhen = r.get('PolyPhen', missing)
	r.EXON = r.get('EXON', missing)
	r.INTRON = r.get('INTRON', missing)
	r.DOMAINS = r.get('DOMAINS', missing)
	r.GMAF = r.get('GMAF', missing)
	r.AFR_MAF = r.get('AFR_MAF', missing)
	r.AMR_MAF = r.get('AMR_MAF', missing)
	r.ASN_MAF = r.get('ASN_MAF', missing)
	r.EAS_MAF = r.get('EAS_MAF', missing)
	r.EUR_MAF = r.get('EUR_MAF', missing)
	r.SAS_MAF = r.get('SAS_MAF', missing)
	r.AA_MAF = r.get('AA_MAF', missing)
	r.EA_MAF = r.get('EA_MAF', missing)
	r.CLIN_SIG = r.get('CLIN_SIG', missing)
	r.SOMATIC = r.get('SOMATIC', missing)
	r.PUBMED = r.get('PUBMED', missing)
	r.MOTIF_NAME = r.get('MOTIF_NAME', missing)
	r.MOTIF_POS = r.get('MOTIF_POS', missing)
	r.HIGH_INF_POS = r.get('HIGH_INF_POS', missing)
	r.MOTIF_SCORE_CHANGE = r.get('MOTIF_SCORE_CHANGE', missing)
	r.IMPACT = r.get('IMPACT', missing)
	r.PICK = r.get('PICK', missing)
	r.VARIANT_CLASS = r.get('VARIANT_CLASS', missing)
	r.TSL = r.get('TSL', missing)
	r.HGVS_OFFSET = r.get('HGVS_OFFSET', missing)
	r.PHENO = r.get('PHENO', missing)
	r.MINIMISED = r.get('MINIMISED', missing)
	r.ExAC_AF = r.get('ExAC_AF', missing)
	r.ExAC_AF_Adj = r.get('ExAC_AF_Adj', missing)
	r.ExAC_AF_AFR = r.get('ExAC_AF_AFR', missing)
	r.ExAC_AF_AMR = r.get('ExAC_AF_AMR', missing)
	r.ExAC_AF_EAS = r.get('ExAC_AF_EAS', missing)
	r.ExAC_AF_FIN = r.get('ExAC_AF_FIN', missing)
	r.ExAC_AF_NFE = r.get('ExAC_AF_NFE', missing)
	r.ExAC_AF_OTH = r.get('ExAC_AF_OTH', missing)
	r.ExAC_AF_SAS = r.get('ExAC_AF_SAS', missing)
	r.GENE_PHENO = r.get('GENE_PHENO', missing)
	r.FILTER = r.get('FILTER', 'PASS')
	r.CONTEXT = r.get('CONTEXT', missing)
	r.src_vcf_id = r.get('src_vcf_id', missing)
	r.tumor_bam_uuid = r.get('tumor_bam_uuid', missing)
	r.normal_bam_uuid = r.get('normal_bam_uuid', missing)
	r.case_id = r.get('case_id', missing)
	r.GDC_FILTER = r.get('GDC_FILTER', missing)
	r.COSMIC = r.get('COSMIC', missing)
	r.MC3_Overlap = r.get('MC3_Overlap', missing)
	r.GDC_Validation_Status = r.get('GDC_Validation_Status', missing)
	r.GDC_Valid_Somatic = r.get('GDC_Valid_Somatic', missing)
	r.vcf_region = r.get('vcf_region', missing)
	r.vcf_info = r.get('vcf_info', missing)
	r.vcf_format = r.get('vcf_format', missing)
	r.vcf_tumor_gt = r.get('vcf_tumor_gt', missing)
	r.vcf_normal_gt = r.get('vcf_normal_gt', missing)

	writer.write(r)
