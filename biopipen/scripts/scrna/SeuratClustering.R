library(Seurat)
library(tidyr)
library(dplyr)

set.seed(8525)

srtfile = {{in.srtobj | quote}}
rdsfile = {{out.rdsfile | quote}}
groupfile = {{out.groupfile | quote}}
envs = {{envs | r}}

sobj = readRDS(srtfile)

sobj = FindNeighbors(sobj)
envs$FindClusters$object = sobj
sobj = do.call(FindClusters, envs$FindClusters)
saveRDS(sobj, file = rdsfile)

groups = NULL
for (ident in unique(Idents(sobj))) {
    df = tibble(cells = WhichCells(sobj, ident = ident)) %>%
        separate("cells", c("Sample", "Cell"), sep="_") %>%
        pivot_wider(
            names_from="Sample",
            values_from="Cell",
            values_fn=function(x) {paste(x, collapse=";")}
        ) %>%
        mutate(Cluster=paste0("Cluster", ident))
    if (is.null(groups)) {
        groups = df
    } else {
        groups = bind_rows(groups, df)
    }
}
groups = groups %>% select(Cluster, everything())

write.table(groups, groupfile, col.names = T, row.names = F, sep = "\t", quote=F)
