source("{{biopipen_dir}}/utils/io.R")
source("{{biopipen_dir}}/utils/plot.R")

# to compile the expressions
library(ComplexHeatmap)

infile = {{in.infile | r}}
annofiles = {{in.annofiles | r}}
outfile = {{out.outfile | r}}
outdir = {{out.outdir | r}}
inopts = {{envs.inopts | r}}
anopts = {{envs.anopts | r}}
drawargs = {{envs.draw | r}}
devpars = {{envs.devpars | r}}
seed = {{envs.seed | r}}

set.seed(seed)

data = read.table.opts(infile, inopts)
annos = lapply(annofiles, function(x) read.table.opts(x, anopts))
if (length(annos) == 1) {
    annos = annos[[1]]
}

# compile the globals
{{envs.globals}}
args = {{envs.args | r}}

hm = plotHeatmap(data, args = args)

do.call(png, c(list(filename=outfile), devpars))
do.call(draw, c(list(hm), drawargs))
dev.off()

saveRDS(hm, file.path(outdir, "heatmap.RDS"))

# export clusters
ro = row_order(hm)
# if no split, treat the whole as one cluster
if (!is.list(ro)) ro = list(ro)
if (is.null(names(ro))) {
    names(ro) = paste0('C', 1:length(ro))
}
rn_orig = rownames(data)
if (is.null(rn_orig)) {
    rn_orig = 1:nrow(data)
}
rclines = c()
for (clname in names(ro)) {
    rclines = c(
        rclines,
        paste("# Cluster:", clname, ', Size:', length(ro[[clname]])),
        paste(rn_orig[ ro[[clname]] ], collapse = ", ")
    )
}
rc_conn = file(file.path(outdir, "row_clusters.txt"))
writeLines(rclines, rc_conn)
close(rc_conn)

co = column_order(hm)
# if no split, treat the whole as one cluster
if (!is.list(co)) co = list(co)
if (is.null(names(co))) {
    names(co) = paste0('C', 1:length(co))
}
cn_orig = colnames(data)
if (is.null(cn_orig)) {
    cn_orig = 1:nrow(data)
}
cclines = c()
for (clname in names(co)) {
    cclines = c(
        cclines,
        paste("# Cluster:", clname, ', Size:', length(co[[clname]])),
        paste(cn_orig[ co[[clname]] ], collapse = ", ")
    )
}
cc_conn = file(file.path(outdir, "col_clusters.txt"))
writeLines(cclines, cc_conn)
close(cc_conn)
