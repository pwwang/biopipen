library(Gviz)
library(GenomicInteractions)

region = unlist(strsplit({{in.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])

irdata = makeGenomicInteractionsFromFile({{in.infile | quote}}, type = {{args.intype | quote}}, experiment_name = {{in.infile | fn | quote}}, description = "Data for {{in.infile | fn}}")

params = list(
	x          = irdata,
	chromosome = chrom,
	name       = {{in.name | quote}}
)
track = do.call(InteractionTrack, params)
displayPars(track) = {{args.params | Rlist}}
saveRDS (track, {{out.outfile | quote}})