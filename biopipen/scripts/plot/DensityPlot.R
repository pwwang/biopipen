library(plotthis)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
val_col <- {{envs.val_col | r}}
devpars <- {{envs.devpars | r}}
args <- {{envs | r}}
args$val_col <- NULL
args$devpars <- NULL

data <- read.table(infile, header=TRUE, sep="\t", check.names=FALSE)
if (is.numeric(val_col)) {
    val_col <- colnames(data)[val_col]
}
val_col <- trimws(val_col)

args$data <- data
args$x <- val_col

p <- do.call(plotthis::DensityPlot, args)

png(outfile, width=devpars$width, height=devpars$height, res=devpars$res)
print(p)
dev.off()
