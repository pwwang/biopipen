source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(immunarch)
library(glue)
library(dplyr)
library(tidyr)
library(tibble)

immfile = {{in.immfile | r}}
sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
metacols = {{envs.metacols | r}}

immdata = readRDS(immfile)
sobj = readRDS(sobjfile)

metadf = do_call(rbind, lapply(seq_len(nrow(immdata$meta)), function(i) {
    # Clones  Proportion   CDR3.aa                       Barcode
    # 5      4 0.008583691 CAVRDTGNTPLVF;CASSEYSNQPQHF   GTTCGGGCACTTACGA-1;TCTCTAAGTACCAGTT-1
    # 6      4 0.008583691 CALTQAAGNKLTF;CASRPEDLRGQPQHF GCTTGAAGTCGGCACT-1;TACTCGCTCCTAAGTG-1
    cldata = immdata$data[[i]][, c(metacols, "Barcode")]
    # # A tibble: 4 × 5
    # Sample                  Patient     Timepoint Tissue
    # <chr>                   <chr>       <chr>     <chr>
    # 1 MC1685Pt011-Baseline-PB MC1685Pt011 Baseline  PB
    mdata = as.list(immdata$meta[i, ])
    for (mname in names(mdata)) {
        assign(mname, mdata[[mname]])
    }

    cldata %>%
        separate_rows(Barcode, sep=";") %>%
        mutate(Barcode = glue("{{envs.prefix}}{Barcode}"))

}))

sobj@meta.data = left_join(
    sobj@meta.data %>% rownames_to_column("Barcode"),
    metadf,
    by = "Barcode"
) %>% column_to_rownames("Barcode")


saveRDS(sobj, outfile)
