#library(Gviz)
#library(GenomicInteractions)

region = unlist(strsplit({{i.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])

infile = {{i.infile | quote}}
intype = {% if args.intype == 'auto' %}{{i.infile | ext | [1:] | quote}}{% else %}{{args.intype | quote}}{% endif %}

if (intype == 'bedpe') {
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
		lines2 = c(lines2, paste0(paste(parts, collapse = '\t'), '\n'))
	}
	rm(lines)
	#infile = textConnection(lines2)
	infile = "{{job.outdir}}/{{i.infile | fn}}.bedpe"
	writeLines(lines2, infile)
	rm(lines2)
} else if (intype == 'bedx') {
	lines = readLines(infile)
	lines2 = NULL
	chr2index    = NA
	start2index  = NA
	end2index    = NA
	name2index   = NA
	strand2index = NA
	for (line in lines) {
		parts = unlist(strsplit(line, '\t', fixed = T))
		if (startsWith(parts[1], "##")) { next }
		if (startsWith(parts[1], "#CHR") || startsWith(parts[1], "CHR")) { # header
			chr2index    = match("CHR2", parts)
			start2index  = match("START2", parts)
			end2index    = match("END2", parts)
			name2index   = match("NAME2", parts)
			strand2index = match("STRAND2", parts)
			next
		} 
		if (is.na(chr2index) || is.na(start2index) || is.na(end2index)) {
			stop('Cannot find CHR2, START2, END2 columns in bedx file')
		}
		name    = if (is.na(name2index)) parts[4] else paste(parts[4], parts[name2index], sep = '-')
		strand2 = if (is.na(strand2index)) '+' else parts[strand2index]
		lines2  = c(lines2, paste0(paste(c(parts[1:3], parts[chr2index], parts[start2index], parts[end2index], name, max(as.character(parts[5]), 1), parts[6], strand2), collapse = '\t'), '\n'))
	}
	rm(lines)
	infile = "{{job.outdir}}/{{i.infile | fn}}.bedpe"
	writeLines(lines2, infile)
	rm(lines2)
	intype = 'bedpe'
}

ret = list()
ret$trackType = 'InteractionTrack'
ret$interactionParams = list(
	fn   = infile,
	type = intype,
	experiment_name = {{i.infile | fn | quote}},
	description = "Interaction data for {{i.infile | fn}}"
)
ret$trackParams = list(
	x          = '',
	chromosome = chrom,
	name       = {{i.name | quote}}
)
ret$trackParams = c(ret$trackParams, {{args.params | Rlist}})

saveRDS (ret, {{o.outfile | quote}})