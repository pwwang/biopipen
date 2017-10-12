if (!exists('txtSampleinfo')) {
	
	txtSampleinfo = function(sfile) {
		mat = read.table (sfile,  header=T, row.names = 1, check.names=F, sep="\t")
		if (length(intersect(c('Patient', 'Group', 'Batch'), colnames(mat))) == 0)
			stop("Headers should be a subset of ['Patient', 'Group', 'Batch']")
		return (mat)
	}
}