# legacy
SampleInfo  = reticulate::import('bioprocs.utils.sampleinfo')$SampleInfo

library(R6)
options(stringsAsFactors = FALSE)

SampleInfo2 = R6Class("SampleInfo2", public = list(

	mat = NULL,
	cnames = c(),

	initialize = function(sifile, checkPaired = FALSE) {
		# 'Sample', 'Patient', 'Group', 'Batch'
		private$read(sifile)
		if (checkPaired && "Patient" %in% self$cnames) {
			for (patient in self$all.patients()) {
				if (length(self$get.samples(by = 'Patient', value = patient)) != 2) {
					stop(sprintf('Expect paired comparisons, but Patient %s has # samples other than 2.', patient))
				}
			}
		}
	},

	select = function(samples) {
		self$mat = self$mat[match(samples, self$mat$Sample),, drop=F]
	},

	as.edger.design = function() {
		private$relevel.group()
		if (self$is.paired()) {
			model.matrix(~ self$mat$Patient + self$mat$Group)
		} else {
			model.matrix(~ self$mat$Group)
		}
	},

	as.deseq2.design = function() {
		private$relevel.group()
		coldata = data.frame(Group = self$mat$Group)
		design  = ~ Group
		if (self$is.paired()) {
			coldata$Patient = self$mat$Patient
			design = ~ Patient + Group
		}
		list(coldata = coldata, design = design)
	},

	as.superfreq.meta = function(datadir = NULL) {
		# BAM, VCF, NAME, INDIVIDUAL, NORMAL, TIMEPOINT
		# Sample  Patient  Group  Batch
		# A.bam   A        Indvd1 NORMAL
		# A.vcf   A        Indvd1 NORMAL
		# B.bam   B        Indvd1 TUMOR
		# B.vcf   B        Indvd1 TUMOR
		ret = data.frame(
			BAM        = character(), VCF    = character(), NAME      = character(),
			INDIVIDUAL = character(), NORMAL = character(), TIMEPOINT = character()
		)
		patients = self$all.patients()
		for (patient in patients) {
			pair = self$get.samples("Patient", patient)
			if (endsWith(pair[1], ".bam")) {
				BAM = pair[1]
				VCF = pair[2]
			} else {
				BAM = pair[2]
				VCF = pair[1]
			}
			BamInfo = self$sample.info(BAM)
			ret = rbind(ret, list(
				BAM        = if (is.null(datadir)) BAM else file.path(datadir, BAM),
				VCF        = if (is.null(datadir)) VCF else file.path(datadir, VCF),
				NAME       = patient,
				INDIVIDUAL = BamInfo$Group,
				NORMAL     = if (toupper(BamInfo$Batch) %in% c("NORMAL", "HEALTHY", "CONTROL")) "YES" else "NO",
				TIMEPOINT  = ""
			))
		}
		ret
	},

	is.paired = function() {
		allpatients = self$all.patients()
		if (is.null(allpatients))
			return (FALSE)
		for (patient in allpatients) {
			if (length(self$get.samples("Patient", patient)) != 2)
				return (FALSE)
		}
		return (TRUE)
	},

	all.patients = function() {
		if (!"Patient" %in% self$cnames)
			return (NULL)
		unique(as.vector(self$mat$Patient))
	},

	all.groups = function() {
		as.vector(levels(relevel(factor(self$mat$Group), self$mat[1, "Group"])))
	},

	get.samples = function(by = NULL, value = NULL) {
		if (!is.null(by) && !by %in% self$cnames) {
			stop(sprintf('%s is not a valid column name.', by))
		}
		if (is.null(by))
			return (as.vector(self$mat$Sample))
		as.vector(self$mat[which(self$mat[, by] == value), "Sample", drop = TRUE])
	},

	sample.info = function(sample, info = NULL) {
		#Get the information of a given sample
		ret = self$mat[which(self$mat$Sample == sample), , drop = FALSE]
		if (is.null(info))
			return (ret)
		return (ret[, info, drop = TRUE])
	}

), private = list(

	read = function(sifile) {
		standard.cnames = c("", "Sample", "Patient", "Group", "Batch")
		self$mat = read.table(sifile, header = TRUE, row.names = NULL, sep = "\t", check.names = FALSE)
		self$cnames = colnames(self$mat)
		if (is.null(self$cnames)) {
			stop('Headers for sample information file is required.')
		}
		if (length(setdiff(self$cnames, standard.cnames)) > 0) {
			stop(sprintf('Headers should be a subset of "%s"', paste(standard.cnames, collapse = ",")))
		}
		if (!"Group" %in% self$cnames) {
			stop("Column 'Group' is required in sample information file.")
		}
		self$cnames = replace(self$cnames, self$cnames == "", "Sample")
		colnames(self$mat) = self$cnames
	},

	relevel.group = function() {
		groups = self$all.groups()
		self$mat$Group = as.factor(self$mat$Group)
		self$mat$Group = relevel(self$mat$Group, groups[2])
	}

))

