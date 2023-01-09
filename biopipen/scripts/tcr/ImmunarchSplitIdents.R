library(stringr)
library(dplyr)
library(immunarch)
library(Seurat)

immdatafile = {{in.immdata | r}}
sobjfile = {{in.sobjfile | r}}
outdir = {{out.outdir | r}}
prefix = {{envs.prefix | r}}

immdata = readRDS(immdatafile)
sobj = readRDS(sobjfile)

get_clusters_for_sample = function(sample) {
    # subset the Seurat object to the sample
    sobj_sample = subset(sobj, subset = {{envs.sample_col}} == sample)
    # get the clusters
    clusters = Idents(sobj_sample)
    # get rid of the prefix
    if (!is.null(prefix) && prefix != "") {
        prefices = str_glue_data(sobj_sample@meta.data, prefix)
        oldnames = names(clusters)
        names(clusters) = sapply(seq_along(clusters), function(i) {
            sub(prefices[i], "", oldnames[i])
        })
    }
    return(clusters)
}


new_data = list()

for (df_i in 1:length(immdata$data)) {
    # The code borrowed from immunarch::select_clusters
    # We need each cluster to be save separately
    source_name <- names(immdata$data)[df_i]
    clusters <- get_clusters_for_sample(source_name)
    # only cells from source name is selected
    df_list <- select_barcodes(immdata$data[[df_i]], clusters, TRUE)

    # Update the data and metadata if everything is OK
    if (length(df_list)) {
        for (new_df_i in 1:length(df_list)) {
            cluster_id <- names(df_list)[new_df_i]
            if (!cluster_id %in% names(new_data)) {
                new_data[[cluster_id]] = list()
            }
            new_data[[cluster_id]][[source_name]] = bind_rows(
                new_data[[cluster_id]][[source_name]],
                df_list[[new_df_i]]
            )
        }
    }
}

for (cluster_id in names(new_data)) {
    meta = immdata$meta
    meta$Cluster = cluster_id
    cldata = list(data = new_data[[cluster_id]], meta = meta)
    saveRDS(cldata, file.path(outdir, paste0(cluster_id, ".rds")))
}
