{% python from os import path %}
library("MatrixEQTL");

{# check if dist is wrongly assign at `args.dist` instead of `args.cisopts.dist`  #}
args = {{args | R}}
if (!is.null(args$dist) && args$cisopts$dist == 0) {
	stop('Set "args.cisopts.dist" instead of "args.dist" to perform cis-eQTL analysis.')
}
snpfile = {{i.snpfile | R}}
expfile = {{i.expfile | R}}
outfile = {{o.outfile | R}}
cisfile = {{o.cisfile | R}}

#snpmatrix = read.table({{i.snpfile | quote}}, sep = "\t", header = T, row.names = 1, check.names = F)
#expmatrix = read.table({{i.expfile | quote}}, sep = "\t", header = T, row.names = 1, check.names = F)
#cnames    = intersect(colnames(snpmatrix), colnames(expmatrix))
#cnames = colnames(expmatrix)
#snpmatrix = snpmatrix[, cnames, drop = F]
#expmatrix = expmatrix[, cnames, drop = F]

snps = SlicedData$new();
snps$fileDelimiter = "\t";       # the TAB character
snps$fileOmitCharacters = "NA";  # denote missing values;
snps$fileSkipRows = 1;           # one row of column labels
snps$fileSkipColumns = 1;        # one column of row labels
snps$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
#snps$CreateFromMatrix( as.matrix(snpmatrix) );
snps$LoadFile( snpfile );

gene = SlicedData$new();
gene$fileDelimiter = "\t";       # the TAB character
gene$fileOmitCharacters = "NA";  # denote missing values;
gene$fileSkipRows = 1;           # one row of column labels
gene$fileSkipColumns = 1;        # one column of row labels
gene$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
#gene$CreateFromMatrix( as.matrix(expmatrix) );
gene$LoadFile( expfile );

cvrt = SlicedData$new();
{% if path.isfile(i.covfile) %}
covmatrix = t(read.table({{i.covfile | quote}}, header = T, row.names = 1, check.names = F))
#covmatrix = covmatrix[, cnames, drop = F]
#cvrt$fileDelimiter = "\t";       # the TAB character
#cvrt$fileOmitCharacters = "NA";  # denote missing values;
#cvrt$fileSkipRows = 1;           # one row of column labels
#cvrt$fileSkipColumns = 1;        # one column of row labels
#cvrt$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
cvrt$CreateFromMatrix( as.matrix(covmatrix) );
{% endif %}

{% if args.cisopts.dist %}
spfile = {{args.cisopts.snppos  | quote}}
gpfile = {{args.cisopts.genepos | quote}}
if (endsWith(spfile, ".bed")) {
	snpspos = read.table(spfile, header = F, sep = '\t', stringsAsFactors = FALSE);
	snpspos = snpspos[, c(4, 1, 2)]
	colnames(snpspos) = c("snp", "chr", "pos")
} else if (endsWith(spfile, ".gff") || endsWith(spfile, ".gtf")) {
	snpspos = read.table(spfile, header = F, sep = '\t', stringsAsFactors = FALSE);
	snpspos = snpspos[, c(9, 1, 4)]
	colnames(snpspos) = c("snp", "chr", "pos")
	snpspos$snp = lapply(snpspos$snp, function(x) for (s in unlist(strsplit(x, '; ', fixed=T))) if (startsWith(s, "gene_id ")) return(substring(s, 9)))
} else {
	snpspos = read.table(spfile, header = TRUE, stringsAsFactors = FALSE);
}
if (endsWith(gpfile, ".bed")) {
	genepos = read.table(gpfile, header = F, sep = '\t', stringsAsFactors = FALSE);
	genepos = genepos[, c(4, 1, 2, 3)]
	colnames(genepos) = c("geneid", "chr", "s1", "s2")
} else if (endsWith(gpfile, ".gff") || endsWith(gpfile, ".gtf")) {
	genepos = read.table(gpfile, header = F, sep = '\t', stringsAsFactors = FALSE);
	genepos = genepos[, c(9, 1, 4, 5)]
	colnames(genepos) = c("geneid", "chr", "s1", "s2")
	genepos$geneid = lapply(genepos$geneid, function(x) for (s in unlist(strsplit(x, '; ', fixed=T))) if (startsWith(s, "gene_id ")) return(substring(s, 9)))
} else {
	genepos = read.table(gpfile, header = TRUE, stringsAsFactors = FALSE);
}

Matrix_eQTL_main(
	snps = snps, 
	gene = gene, 
	cvrt = cvrt, 
	output_file_name = outfile, 
	pvOutputThreshold = {{args.pval}},
	useModel = {{args.model | lambda x: 'r:' + x | R}}, 
	errorCovariance = numeric(), 
	verbose=T, 
	output_file_name.cis  = cisfile,
	pvOutputThreshold.cis = {{args.cisopts.cispv}},
	snpspos = snpspos, 
	genepos = genepos,
	cisDist = {{args.cisopts.dist}},
	noFDRsaveMemory = !{{args.fdr | R}}
)
# if trans-eqtl analysis is not done (args.pval = 0)
if (!file.exists(outfile))
	file.create(outfile)
{% else %}
Matrix_eQTL_engine(
	snps = snps, 
	gene = gene, 
	cvrt = cvrt, 
	output_file_name = outfile, 
	pvOutputThreshold = {{args.pval}}, 
	useModel = {{args.model | lambda x: 'r:' + x | R}}, 
	errorCovariance = numeric(), 
	verbose = T,
	noFDRsaveMemory = !{{args.fdr | R}}
)
file.create(cisfile)
{% endif %}