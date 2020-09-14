library(philentropy)
{{rimport}}('__init__.r')

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R}}
method  = {{args.method | R}}
transfm = {{args.transfm | R}}
na      = {{args.na | R}}

indata = read.table.inopts(infile, inopts)
indata[is.na(indata)] = na

d = distance(indata, method = method, test.na = F)
if (transfm != FALSE) {
	ma = max(d)
	mi = min(d)
	scaled = 1 - (d-mi) / (ma-mi)
	if (transfm == 'scale') {
		d = scaled
	} else if (transfm == 'similarity' || transfm == 'sim') {
		d = 1 - scaled
	}
}
rownames(d) = rownames(indata)
colnames(d) = rownames(indata)
write.table(round(d, 3), outfile, col.names = TRUE, row.names = TRUE, sep = "\t", quote = F)
