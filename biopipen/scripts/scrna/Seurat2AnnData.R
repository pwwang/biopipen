library(biopipen.utils)

sobjfile <- {{in.sobjfile | quote}}
outfile <- {{out.outfile | quote}}
assay <- {{envs.assay | r}}

ConvertSeuratToAnnData(sobjfile, outfile = outfile, assay = assay)
