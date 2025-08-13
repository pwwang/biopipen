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

log <- get_logger()

log$info("Loading scRepertoire object ...")
screp <- read_obj(screpfile)

log$info("Loading Seurat object ...")
srtobj <- read_obj(srtobjfile)

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
