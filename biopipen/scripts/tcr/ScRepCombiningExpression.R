library(scRepertoire)
library(Seurat)
library(biopipen.utils)

screpfile <- {{in.screpfile | r}}
srtobjfile <- {{in.srtobj | r}}
outfile <- {{out.outfile | r}}
cloneCall <- {{envs.cloneCall | r}}
chain <- {{envs.chain | r}}
group_by <- {{envs.group_by | default: envs["group-by"] | default: None | r}}
proportion <- {{envs.proportion | r}}
filterNA <- {{envs.filterNA | r}}
cloneSize <- {{envs.cloneSize | r}}
addLabel <- {{envs.addLabel | r}}
imm_cell_id_trans <- {{envs.imm_cell_id_trans | r}}
rna_cell_id_trans <- {{envs.rna_cell_id_trans | r}}

if (!is.null(imm_cell_id_trans)) {
    if (!is.character(imm_cell_id_trans) && !is.function(imm_cell_id_trans)) {
        stop("`envs.imm_cell_id_trans` must be a string of R function")
    }
    if (is.character(imm_cell_id_trans)) {
        imm_cell_id_trans <- eval(parse(text = imm_cell_id_trans))
    }
}

if (!is.null(rna_cell_id_trans)) {
    if (!is.character(rna_cell_id_trans) && !is.function(rna_cell_id_trans)) {
        stop("`envs.rna_cell_id_trans` must be a string of R function")
    }
    if (is.character(rna_cell_id_trans)) {
        rna_cell_id_trans <- eval(parse(text = rna_cell_id_trans))
    }
}

log <- get_logger()

log$info("Loading scRepertoire object ...")
screp <- read_obj(screpfile)
if (!is.null(imm_cell_id_trans)) {
    log$info("Transforming cell IDs in scRepertoire object ...")
    screp <- lapply(screp, function(x) {
        x$barcode <- imm_cell_id_trans(x$barcode)
        return(x)
    })
}

log$info("Loading Seurat object ...")
srtobj <- read_obj(srtobjfile)
if (!is.null(rna_cell_id_trans)) {
    log$info("Transforming cell IDs in Seurat object ...")
    srtobj <- RenameCells(srtobj, new.names = rna_cell_id_trans(Cells(srtobj)))
}

log$info("Combining expression data ...")

obj <- combineExpression(
    input.data = screp,
    sc.data = srtobj,
    cloneCall = cloneCall,
    chain = chain,
    group.by = group_by,
    proportion = proportion,
    filterNA = filterNA,
    cloneSize = unlist(cloneSize),
    addLabel = addLabel
)
obj$VDJ_Presence <- !is.na(obj$CTaa)

log$info("Saving combined object ...")
save_obj(obj, outfile)
