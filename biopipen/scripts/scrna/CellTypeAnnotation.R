library(Seurat)
library(biopipen.utils)
set.seed(8525)

backup_col <- {{envs.backup_col | r}}


merge_clusters_with_same_labels <- function(sobj, newcol = NULL) {
    if (is.null(newcol)) {
        newcol <- biopipen.utils::GetIdentityColumn(sobj)
        sobj@meta.data[[newcol]] <- sub("\\.\\d+$", "", sobj@meta.data[[newcol]])
        Idents(sobj) <- newcol
    } else {
        sobj@meta.data[[newcol]] <- sub("\\.\\d+$", "", sobj@meta.data[[newcol]])
    }

    sobj
}

rename_idents <- function(sobj, ident_col, mapping) {
    orig_ident_col <- biopipen.utils::GetIdentityColumn(sobj)
    if (!identical(ident_col, orig_ident_col)) {
        Idents(sobj) <- ident_col
        mapping$object <- sobj
        sobj <- do_call(RenameIdents, mapping)
    } else {
        if (!is.null(backup_col)) {
            sobj@meta.data[[backup_col]] <- Idents(sobj)
        }
        mapping$object <- sobj
        sobj <- do_call(RenameIdents, mapping)
    }
    Idents(sobj) <- ident_col
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
