library("MatrixEQTL");

snps = SlicedData$new();
snps$fileDelimiter = "\t";       # the TAB character
snps$fileOmitCharacters = "NA";  # denote missing values;
snps$fileSkipRows = 1;           # one row of column labels
snps$fileSkipColumns = 1;        # one column of row labels
snps$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps$LoadFile( {{in.snpfile | quote}} );

gene = SlicedData$new();
gene$fileDelimiter = "\t";       # the TAB character
gene$fileOmitCharacters = "NA";  # denote missing values;
gene$fileSkipRows = 1;           # one row of column labels
gene$fileSkipColumns = 1;        # one column of row labels
gene$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
gene$LoadFile( {{in.expfile | quote}} );

cvrt = SlicedData$new();
{% if in.covfile %}
cvrt$fileDelimiter = "\t";       # the TAB character
cvrt$fileOmitCharacters = "NA";  # denote missing values;
cvrt$fileSkipRows = 1;           # one row of column labels
cvrt$fileSkipColumns = 1;        # one column of row labels
cvrt$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
cvrt$LoadFile( {{in.covfile | quote}} );
{% endif %}

{% if args.cisopts.dist %}
snpspos = read.table({{args.cisopts.snppos  | quote}}, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table({{args.cisopts.genepos | quote}}, header = TRUE, stringsAsFactors = FALSE);
Matrix_eQTL_main(
	snps = snps, 
	gene = gene, 
	cvrt = cvrt, 
	output_file_name = {{out.outfile | quote}}, 
	pvOutputThreshold = {{args.pval}},
	useModel = {{args.model | lambda x: 'r:' + x | R}}, 
	errorCovariance = numeric(), 
	verbose=T, 
	output_file_name.cis  = {{out.cisfile | quote}},
	pvOutputThreshold.cis = {{args.cisopts.cispv}},
	snpspos = snpspos, 
	genepos = genepos,
	cisDist = {{args.cisopts.dist}},
	noFDRsaveMemory = !{{args.fdr | R}}
)
{% else %}
Matrix_eQTL_engine(
	snps = snps, 
	gene = gene, 
	cvrt = cvrt, 
	output_file_name = {{out.outfile | quote}}, 
	pvOutputThreshold = {{args.pval}}, 
	useModel = {{args.model | lambda x: 'r:' + x | R}}, 
	errorCovariance = numeric(), 
	verbose = T,
	noFDRsaveMemory = !{{args.fdr | R}}
)
file.create({{out.cisfile | quote}})
{% endif %}