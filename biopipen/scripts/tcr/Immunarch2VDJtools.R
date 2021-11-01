library(immunarch)
library(dplyr)
library(tidyr)
library(stringr)

immfile = {{in.immdata | quote}}
outdir = {{out.outdir | quote}}

immdata = readRDS(immfile)

for (sample in names(immdata$data)) {
    # see https://vdjtools-doc.readthedocs.io/en/master/input.html
    # for the input format

    df = immdata$data[[sample]] %>%
        transmute(
            contig_id=ContigID,
            count=Clones,
            frequency=Proportion,
            CDR3nt=CDR3.nt,
            CDR3aa=CDR3.aa,
            V=V.name,
            D=D.name,
            J=J.name
        ) %>%
        tibble::as_tibble() %>%
        rowwise() %>%
        mutate(
            contig_id=paste(
                unlist(strsplit(contig_id, ";"))[1:(str_count(CDR3nt, ";")+1)],
                collapse=";"
            )
        ) %>%
        separate_rows(contig_id, CDR3nt, CDR3aa, V, D, J, sep=";") %>%
        dplyr::select(-1)

    outfile = file.path(outdir, paste0(sample, ".vdjtools.txt"))
    write.table(
        df, outfile,
        sep="\t", quote=F, row.names = F, col.names=T
    )
}
