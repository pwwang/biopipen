{{"__init__.R" | rimport}}
options(stringsAsFactors=FALSE)
library("MatrixEQTL")

snpfile = {{i.snpfile | R}}
expfile = {{i.expfile | R}}
covfile = {{i.covfile | R}}
outfile = {{o.outfile | R}}
alleqtl = {{o.alleqtl | R}}

model = {{args.model | R}}
pval = {{args.pval | R}}
transp = {{args.transp | ?isinstance: bool | =_ | ?=: 1e-5 | $R}}
fdr = {{args.fdr | R}}
snppos = {{args.snppos | R}}
genepos = {{args.genepos | R}}
dist = {{args.dist | R}}

if (!model %in% c("modelANOVA", "modelLINEAR", "modelLINEAR_CROSS")) {
	stop("Model should be one of modelANOVA, modelLINEAR, or modelLINEAR_CROSS")
}
model = get(model)

trans_enabled = is.true(transp)
cis_enabled = is.true(snppos) && is.true(genepos) && dist > 0

# if trans is disabled, all files needed for cis should be provided
if (!trans_enabled && !cis_enabled) {
	stop(paste("Trans-eQTL is disabled but parameter (snppos, genepos or dist)",
		 "needed for cis-eQTL is not provided."))
}

snps = SlicedData$new();
snps$fileDelimiter = "\t";       # the TAB character
snps$fileOmitCharacters = "NA";  # denote missing values;
snps$fileSkipRows = 1;           # one row of column labels
snps$fileSkipColumns = 1;        # one column of row labels
snps$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps$LoadFile( snpfile );

gene = SlicedData$new();
gene$fileDelimiter = "\t";       # the TAB character
gene$fileOmitCharacters = "NA";  # denote missing values;
gene$fileSkipRows = 1;           # one row of column labels
gene$fileSkipColumns = 1;        # one column of row labels
gene$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
gene$LoadFile( expfile );

cvrt = SlicedData$new();
if (is.true(covfile) && file.exists(covfile)) {
	covmatrix = t(read.table.inopts(covfile, list(cnames=TRUE, rnames=TRUE)))
	cvrt$CreateFromMatrix( as.matrix(covmatrix) )
}

engine_params = list()
engine_params$snps = snps
engine_params$gene = gene
engine_params$cvrt = cvrt
engine_params$output_file_name = ifelse(trans_enabled, alleqtl, NULL)
engine_params$pvOutputThreshold = ifelse(trans_enabled, transp, 0)
engine_params$useModel = model
engine_params$errorCovariance = numeric()
engine_params$verbose = TRUE
engine_params$noFDRsaveMemory = !fdr

noq = function(s) {
	gsub('^\"|\"$', "", s)
}

if (cis_enabled) {
	if (endsWith(snppos, ".bed")) {
		snppos_data = read.table.inopts(snppos,
										list(cnames=FALSE, rnames=FALSE))
		snppos_data = snppos_data[, c(4, 1, 2)]
		colnames(snppos_data) = c("snp", "chr", "pos")
	} else if (endsWith(snppos, ".gff") || endsWith(snppos, ".gtf")) {
		snppos_data = read.table.inopts(snppos,
										list(cnames=FALSE, rnames=FALSE));
		snppos_data = snppos_data[, c(9, 1, 4)]
		colnames(snppos_data) = c("snp", "chr", "pos")
		snppos_data$snp = unlist(lapply(snppos_data$snp, function(x) {
			for (s in unlist(strsplit(x, '; ', fixed=T))) {
				if (startsWith(s, "snp_id "))
					return(noq(substring(s, 8)))
				else if (startsWith(s, "rs_id "))
					return(noq(substring(s, 7)))
				else if (startsWith(s, "rs "))
					return(noq(substring(s, 4)))
			}
		}))
	} else if (endsWith(snppos, ".vcf") || endsWith(snppos, ".vcf.gz")) {
		snppos_data = read.table.inopts(snppos,
										list(cnames=FALSE, rnames=FALSE))
		snppos_data = snppos_data[, c(3, 1, 2)]
		colnames(snppos_data) = c("snp", "chr", "pos")
	} else {
		snppos_data = read.table.inopts(snppos, list(cnames=TRUE))
		colnames(snppos_data) = c("snp", "chr", "pos")
	}

	if (endsWith(genepos, ".bed")) {
		genepos_data = read.table.inopts(genepos,
										 list(cnames=FALSE, rnames=FALSE))
		genepos_data = genepos_data[, c(4, 1:3)]
		colnames(genepos_data) = c("geneid", "chr", "s1", "s2")
	} else if (endsWith(genepos, ".gff") || endsWith(genepos, ".gtf")) {
		genepos_data = read.table.inopts(genepos,
										 list(cnames=FALSE, rnames=FALSE))
		genepos_data = genepos_data[, c(9, 1, 4, 5)]
		colnames(genepos_data) = c("geneid", "chr", "s1", "s2")
		genepos_data$geneid = noquote(unlist(lapply(genepos_data$geneid, function(x) {
			for (s in unlist(strsplit(x, '; ', fixed=T))) {
				if (startsWith(s, "gene_id "))
					return(noq(substring(s, 9)))
			}
		})))
	} else {
		genepos_data = read.table(genepos, header = TRUE, stringsAsFactors = FALSE);
		colnames(genepos_data) = c("geneid", "chr", "s1", "s2")
	}

	engine_params$output_file_name.cis = outfile
	engine_params$pvOutputThreshold.cis = pval
	engine_params$cisDist = dist
	engine_params$snpspos = snppos_data
	engine_params$genepos = genepos_data
	do.call(Matrix_eQTL_main, engine_params)
} else {
	do.call(Matrix_eQTL_engine, engine_params)
}

# if trans-eqtl analysis is not done (args.pval = 0)
if (!file.exists(outfile))
	file.create(outfile)
if (!file.exists(alleqtl))
	file.create(alleqtl)
