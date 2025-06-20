library(scRepertoire)
library(Seurat)
library(biopipen.utils)

screpfile <- {{in.screpfile | r}}
srtobjfile <- {{in.srtobj | r}}
outfile <- {{out.outfile | r}}
cloneCall <- {{envs.cloneCall | r}}
chain <- {{envs.chain | r}}
group.by <- {{envs["group-by"] | r}}
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
    group.by = group.by,
    proportion = proportion,
    filterNA = filterNA,
    cloneSize = unlist(cloneSize),
    addLabel = addLabel
)

log$info("Saving combined object ...")
save_obj(obj, outfile)
