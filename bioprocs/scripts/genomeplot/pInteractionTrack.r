library(Gviz)
library(GenomicInteractions)

region = unlist(strsplit({{in.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])

infile = {{in.infile | quote}}
{% if args.intype | lambda x: x == 'bedpe' %}
# check bedpe format
# normally it should be: (see http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
# chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2[	other fields]
# chr1	3427000	2343059	chr1	2351020	2351159	42401_42410	3.783	+	-	143	79	3
# chr16	3043880	3044279	chr16	3056380	3056559	12898_12909	2.344	-	*	492	89	2
# !!! No header
# !!! Strand must be +, - or *

# Try to remove header and convert strand . to * (this format is usually from Encode)
lines = readLines(infile)
# has header
if (startsWith(lines[1], "chrom")) {
	lines = lines[-1]
}
lines2 = NULL
for (line in lines) {
	parts = unlist(strsplit(line, '\t', fixed = T))
	# filter out interaction for other chromosomes, cuz you can't plot it!
	# TODO: Also do for other input file types! Dramatically save time!
	if (parts[1] != chrom) { next } 
	if (parts[9] == '.')  parts[9]  = '*'
	if (parts[10] == '.') parts[10] = '*'
	#lines[i] = paste(parts, collapse = '\t')
	lines2 = c(lines2, paste(parts, collapse = '\t'))
}
rm(lines)
infile = textConnection(lines2)
{% endif %}

irdata = makeGenomicInteractionsFromFile(infile, type = {{args.intype | quote}}, experiment_name = {{in.infile | fn | quote}}, description = "Data for {{in.infile | fn}}")

params = list(
	x          = irdata,
	chromosome = chrom,
	name       = {{in.name | quote}}
)
track = do.call(InteractionTrack, params)
displayPars(track) = {{args.params | Rlist}}
saveRDS (track, {{out.outfile | quote}})