set.seed(8525)

merge_clusters_with_same_labels <- function(sobj, newcol) {
    if (is.null(newcol)) {
        newcol <- biopipen.utils::GetIdentityColumn(sobj)
        sobj@meta.data[[newcol]] <- sub("\\.\\d+$", "", sobj@meta.data[[newcol]])
        Idents(sobj) <- newcol
    } else {
        sobj@meta.data[[newcol]] <- sub("\\.\\d+$", "", sobj@meta.data[[newcol]])
    }

    sobj
}

{% if envs.tool == "hitype" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-hitype.R" %}
{% elif envs.tool == "sctype" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-sctype.R" %}
{% elif envs.tool == "sccatch" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-sccatch.R" %}
{% elif envs.tool == "celltypist" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-celltypist.R" %}
{% elif envs.tool == "direct" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-direct.R" %}
{% else %}
stop("Unknown tool: {{envs.tool}}")
{% endif %}
