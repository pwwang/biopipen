library(immunarch)
library(dplyr)
library(tidyr)
library(stringr)

immfile = {{in.immdata | r}}
outdir = {{out.outdir | r}}

immdata = readRDS(immfile)

for (sample in names(immdata$data)) {
    # see https://vdjtools-doc.readthedocs.io/en/master/input.html
    # for the input format

    df = immdata$data[[sample]] %>%
        transmute(
            count=Clones,
            frequency=Proportion,
            CDR3nt=CDR3.nt,
            CDR3aa=CDR3.aa,
            V=V.name,
            D=D.name,
            J=J.name
        )

    outfile = file.path(outdir, paste0(sample, ".vdjtools.txt"))
    write.table(
        df, outfile,
        sep="\t", quote=F, row.names = F, col.names=T
    )
}
