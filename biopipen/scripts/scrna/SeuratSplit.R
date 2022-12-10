library(Seurat)

srtobjfile = {{in.srtobj | r}}
inby = {{in.by | r}}
outdir = {{out.outdir | r}}
by = {{envs.by | r}}
recell = NULL
{% if envs.recell %}
recell = {{envs.recell}}
{% endif %}

if (is.null(inby) && is.null(by)) {
    stop("Either `in.by` or `envs.by` must be specified")
}

if (!is.null(inby)) {
    by = inby
}

srtobj = readRDS(srtobjfile)
objlist = SplitObject(srtobj, split.by = by)

for (i in 1:length(objlist)) {
    obj = objlist[[i]]
    name = obj@meta.data[[by]][1]

    if (!is.null(recell)) {
        oldids = rownames(obj@meta.data)
        bydata = obj@meta.data[[by]]
        newids = recell(oldids, bydata)
        obj = RenameCells(obj, new.names = newids)
    }
    saveRDS(obj, file.path(outdir, paste0(name, ".rds")))
}
