{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "single_cell.R" | source_r }}

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
assay <- {{envs.assay | r}}

seurat_to_anndata(sobjfile, outfile, assay, log_info)
