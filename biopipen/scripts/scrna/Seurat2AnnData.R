library(biopipen.utils)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
assay <- {{envs.assay | r}}

ConvertSeuratToAnnData(sobjfile, outfile = outfile, assay = assay)
