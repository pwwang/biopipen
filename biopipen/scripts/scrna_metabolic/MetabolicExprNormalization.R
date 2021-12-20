source("{{biopipen_dir}}/utils/rnaseq.R")

library(scran)

sceobjfile = {{in.sceobj | r}}
config = {{in.configfile | read | toml_loads | r}}
outfile = {{out.outfile | r}}
dropout_cutoff = {{envs.dropout | r}}
refexon = {{envs.refexon | r}}

selected_impute_sce = readRDS(sceobjfile)
cell_types = unique(selected_impute_sce[[config$grouping_name]])

gene_select_mat = matrix(
    FALSE,
    nrow=nrow(selected_impute_sce),
    ncol=length(cell_types),
    dimnames = list(rownames(selected_impute_sce), cell_types)
)

for(c in cell_types){
  each_sce = selected_impute_sce[,selected_impute_sce[[config$grouping_name]] == c]
  each_exp = assay(each_sce,"exprs")
  dropout_rate = apply(each_exp, 1, function(x) sum(x>0)/ncol(each_exp))
  select = dropout_rate >= dropout_cutoff
  gene_select_mat[select, c] = TRUE
}

print("The number of genes selected:")
print(sum(rowSums(gene_select_mat) >= length(cell_types)))

low_dropout_genes = rownames(gene_select_mat)[
    rowSums(gene_select_mat) >= length(cell_types)
]

selected_impute_tpm = tpm(selected_impute_sce)
genelen = glenFromGFFExons(refexon)
genelen = as.numeric(as.vector(unlist(genelen[rownames(selected_impute_sce)])))
selected_impute_counts = sweep(selected_impute_tpm, 1, genelen, FUN = "*")

scran_sf = computeSumFactors(
    SingleCellExperiment(list(counts=selected_impute_counts[low_dropout_genes,])),
    clusters=selected_impute_sce[[config$grouping_name]]
)
summary(scran_sf)
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / scran_sf$sizeFactor)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm+1)

sce_norm = SingleCellExperiment(
    assays = list(tpm=selected_impute_tpm_norm, exprs=selected_impute_exp_norm),
    colData = colData(selected_impute_sce),
    rowData = rowData(selected_impute_sce)
)

saveRDS(sce_norm, file = outfile)
