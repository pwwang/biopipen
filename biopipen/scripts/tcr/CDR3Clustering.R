library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(biopipen.utils)

screpfile <- {{in.screpfile | r}}
outdir <- normalizePath({{job.outdir | r}})
outfile <- {{out.outfile | r}}

tool <- {{envs.tool | r}}
python <- {{envs.python | r}}
within_sample <- {{envs.within_sample | r}}
args <- {{envs.args | r}}
chain <- {{envs.chain | r}}
type <- {{envs.type | r}}

setwd(outdir)

log <- get_logger()

log$info("Reading input file ...")
obj <- read_obj(screpfile)
is_seurat <- inherits(obj, "Seurat")


get_type <- function() {
    if (!is_seurat) {
        for (sample in names(obj)) {
            for (gene in obj[[sample]]$CTgene) {
                if (grepl("^TRB", gene) || grepl("^TRG", gene) || grepl("^TRA", gene) || grepl("^TRD", gene)) {
                    return("TCR")
                } else if (grepl("^IGH", gene) || grepl("^IGK", gene) || grepl("^IGL", gene)) {
                    return("BCR")
                }
            }
        }
    } else {
        for (gene in obj@meta.data$CTgene) {
            if (grepl("^TRB", gene) || grepl("^TRG", gene) || grepl("^TRA", gene) || grepl("^TRD", gene)) {
                return("TCR")
            } else if (grepl("^IGH", gene) || grepl("^IGK", gene) || grepl("^IGL", gene)) {
                return("BCR")
            }
        }
    }
    stop("Cannot determine the type of the data (TCR or BCR). Please set envs.type to 'TCR' or 'BCR'.")
}

if (type == "auto") {
    type <- get_type()
    log$info("Auto-detected data type: {type}")
}

get_cdr3aa_df <- function() {
    if (!is_seurat) {
        out <- NULL
        for (sample in names(obj)) {
            df <- data.frame(
                Sample = sample,
                Barcode = obj[[sample]]$barcode
            )
            if (chain == "both") {
                df$CDR3.aa <- obj[[sample]]$CTaa
            } else if ((type == "BCR" && chain == "heavy") || (type == "TCR" && chain == "light")) {
                df$CDR3.aa <- obj[[sample]]$cdr3_aa1
            } else if ((type == "BCR" && chain == "light") || (type == "TCR" && chain == "heavy")) {
                df$CDR3.aa <- obj[[sample]]$cdr3_aa2
            } else {
                stop(paste("Unknown chain:", chain, "for", type))
            }
            out <- rbind(out, df)
        }
    } else {
        out <- obj@meta.data
        out$Barcode <- rownames(out)
        out <- out %>% filter(!is.na(CTaa))
        if (grepl("_", out$CTaa[1])) {
            if (chain == "both") {
                out$CDR3.aa <- out$CTaa
            } else {
                out <- separate(out, CTaa, into = c("first", "second"), sep = "_")
                if ((type == "BCR" && chain == "heavy") || (type == "TCR" && chain == "light")) {
                    out$CDR3.aa <- out$first
                } else if ((type == "BCR" && chain == "light") || (type == "TCR" && chain == "heavy")) {
                    out$CDR3.aa <- out$second
                } else {
                    stop(paste("Unknown chain:", chain, "for", type))
                }
            }
        } else {
            out$CDR3.aa <- out$CTaa
        }
        out <- select(out, Sample, Barcode, CDR3.aa)
    }

    # Sample, Barcode, CDR3.aa
    out
}
cdr3aa_df = get_cdr3aa_df()

prepare_clustcr = function(clustcr_dir) {
    clustering_args = ""
    for (name in names(args)) {
        value = args[[name]]
        if (is.logical(value)) {
            value = tools::toTitleCase(as.character(value))
        } else if (is.character(value)) {
            value = paste0("'", value, "'")
        }
        clustering_args = paste(name, "=", value)
    }
    clustcr_source = '
import sys
import atexit

import pandas as pd
from scipy import sparse as scipy_sparse


@atexit.register
def clustcr_exit():
    import pandas as pd
    import numpy
    import scipy
    import sklearn
    import matplotlib
    sys.stderr.write("Session info:\\n")
    sys.stderr.write(f"- pandas: {pd.__version__}\\n")
    sys.stderr.write(f"- numpy: {numpy.__version__}\\n")
    sys.stderr.write(f"- scipy: {scipy.__version__}\\n")
    sys.stderr.write(f"- sklearn: {sklearn.__version__}\\n")
    sys.stderr.write(f"- matplotlib: {matplotlib.__version__}\\n")


# Monkey-patch scipy.sparse.isspmatrix to adopt latest scipy v1.14
# If not, an error is raised:
#   numpy.linalg.LinAlgError: 0-dimensional array given.
#   Array must be at least two-dimensional
scipy_sparse.isspmatrix = lambda x: isinstance(
    x,
    (
        scipy_sparse.spmatrix,
        scipy_sparse.csr_array,
        scipy_sparse.csr_matrix,
        scipy_sparse.csc_array,
        scipy_sparse.csc_matrix,
    ),
)


import clustcr  # noqa: #402

clustcr_dir, clustcr_infile = sys.argv[1:3]
cdr3df = pd.read_csv(clustcr_infile, index_col=None)
cdr3 = cdr3df.iloc[:, 0]

clustering = clustcr.Clustering()
output = clustering.fit(cdr3)
output.clusters_df.to_csv(clustcr_dir + "/clusters.txt", sep="\t", index=False)
'
    clustcr_file = file.path(clustcr_dir, "_clustcr.py")
    cat(sprintf(clustcr_source, clustering_args), file=clustcr_file)
    clustcr_file
}

clean_clustcr_output = function(clustcr_outfile) {
    clustcr_out = read.delim2(clustcr_outfile, header=TRUE, row.names = NULL)
    colnames(clustcr_out) = c("CDR3.aa", "CDR3_Cluster")
    out = left_join(cdr3aa_df, distinct(clustcr_out), by=c(cdr3seq4clustering = "CDR3.aa")) %>%
        mutate(
            CDR3_Cluster = if_else(
                is.na(CDR3_Cluster),
                paste0("S_", row_number()),
                paste0("M_", as.character(CDR3_Cluster))
            )
        )

    if (within_sample) {
        out <- mutate(out, CDR3_Cluster = paste0(Sample, ".", CDR3_Cluster))
    }

    # This join would result in more rows than dplyr can handle
    # left_join(cdr3aa_df, out, by = "CDR3.aa")
    out <- out[match(cdr3aa_df$CDR3.aa, out$CDR3.aa), , drop=FALSE]
    cbind(cdr3aa_df, out[, setdiff(colnames(out), "CDR3.aa"), drop=FALSE])
}

run_clustcr = function() {
    log$info("Running ClusTCR ...")
    clustcr_dir = file.path(outdir, "ClusTCR_Output")
    dir.create(clustcr_dir, showWarnings = FALSE)
    clustcr_file = prepare_clustcr(clustcr_dir)
    clustcr_input = prepare_input()
    clustcr_cmd = paste(
        python,
        clustcr_file,
        clustcr_dir,
        clustcr_input
    )
    print("Running:")
    print(clustcr_cmd)
    log$debug("- Running command: {clustcr_cmd}")
    rc = system(clustcr_cmd)
    if (rc != 0) {
        quit(status=rc)
    }
    clustcr_outfile = file.path(clustcr_dir, "clusters.txt")
    clean_clustcr_output(clustcr_outfile)
}

prepare_giana = function() {
    biopipen_dir <- get_biopipen_dir(python)
    giana_srcdir = file.path(biopipen_dir, "scripts", "tcr", "GIANA")

    # # The source code of GIANA is downloaded now to giana_srcdir
    # giana_file = file.path(giana_srcdir, "GIANA.py")
    # giana4_file = file.path(giana_srcdir, "GIANA4.py")
    # giana_query = file.path(giana_srcdir, "query.py")
    # giana_trbv = file.path(giana_srcdir, "Imgt_Human_TRBV.fasta")
    # if (!file.exists(giana_file)) {
    #     download.file(paste(giana_repo, "GIANA4.1.py", sep="/"), giana_file)
    #     download.file(paste(giana_repo, "GIANA4.py", sep="/"), giana4_file)
    #     download.file(paste(giana_repo, "query.py", sep="/"), giana_query)
    #     download.file(paste(giana_repo, "Imgt_Human_TRBV.fasta", sep="/"), giana_trbv)
    # }

    giana_srcdir
}

prepare_input = function() {
    cdr3aa_df$cdr3seq4clustering <<- gsub("[^A-Z]", "", cdr3aa_df$CDR3.aa)  # Remove non-amino acid characters
    cdr3 <- unique(cdr3aa_df$cdr3seq4clustering)

    # cdr3 = distinct(cdr3, aminoAcid, vMaxResolved)

    cdr3file = file.path(outdir, "cdr3.csv")
    write.table(
        data.frame(CDR3.aa=cdr3),
        cdr3file,
        row.names=FALSE, col.names=TRUE, quote=FALSE
    )
    cdr3file
}

clean_giana_output = function(giana_outfile) {
    # generate an output file with columns:
    # CDR3.aa, CDR3_Cluster, V.name, Sample
    # If sequence doesn't exist in the input file,
    # Then a unique cluster id is assigned to it.
    giana_out = read.delim2(giana_outfile, header=FALSE, comment.char = "#", row.names = NULL)[, 1:2, drop=FALSE]
    colnames(giana_out) = c("CDR3.aa", "CDR3_Cluster")
    out = left_join(cdr3aa_df, distinct(giana_out), by=c(cdr3seq4clustering = "CDR3.aa")) %>%
        mutate(
            CDR3_Cluster = if_else(
                is.na(CDR3_Cluster),
                paste0("S_", row_number()),
                paste0("M_", as.character(CDR3_Cluster))
            )
        )

    if (within_sample) {
        out <- mutate(out, CDR3_Cluster = paste0(Sample, ".", CDR3_Cluster))
    }

    # This join would result in more rows than dplyr can handle
    # left_join(cdr3aa_df, out, by = "CDR3.aa")
    out <- out[match(cdr3aa_df$CDR3.aa, out$CDR3.aa), , drop=FALSE]
    cbind(cdr3aa_df, out[, setdiff(colnames(out), "CDR3.aa"), drop=FALSE])
}

run_giana = function() {
    log$info("Running GIANA ...")
    giana_srcdir = prepare_giana()
    giana_input = prepare_input()
    giana_outdir = file.path(outdir, "GIANA_Output")
    dir.create(giana_outdir, showWarnings = FALSE)
    args_str = ""
    for (argname in names(args)) {
        argvalue = args[[argname]]
        if (!startsWith(argname, "-")) {
            if (nchar(argname) == 1) {
                argname = paste0("-", argname)
            } else {
                argname = paste0("--", argname)
            }
        }
        if (isTRUE(argvalue) || toupper(as.character(argvalue)) == "TRUE") {
            argvalue = ""
        } else {
            argvalue = as.character(argvalue)
        }
        args_str = paste(args_str, argname, argvalue)
    }
    giana_cmd = paste(
        python,
        file.path(giana_srcdir, "GIANA.py"),
        "-f", giana_input,
        "-o", giana_outdir,
        "-v", # TRBV mutation not supported
        args_str
    )
    print("Running:")
    print(giana_cmd)
    log$debug("- Running command: {giana_cmd}")
    rc = system(giana_cmd)
    if (rc != 0) {
        quit(status=rc)
    }
    giana_outfile = file.path(giana_outdir, "cdr3--RotationEncodingBL62.txt")
    clean_giana_output(giana_outfile)
}

attach_to_obj = function(obj, out) {
    out <- as.data.frame(out)
    rownames(out) <- out$Barcode
    if (is_seurat) {
        # Attach results to Seurat object
        obj@meta.data$CDR3_Cluster <- out[rownames(obj@meta.data), "CDR3_Cluster"]
    } else {
        # Attach results to the list of data frames
        for (sample in names(obj)) {
            sout <- filter(out, Sample == sample)
            obj[[sample]]$CDR3_Cluster <- sout[obj[[sample]]$barcode, "CDR3_Cluster"]
        }
    }
    obj
}


if (tolower(tool) == "clustcr") {
    out = run_clustcr()
} else if (tolower(tool) == "giana") {
    out = run_giana()
} else {
    stop(paste("Unknown tool:", tool))
}

log$info("Attaching results to the input object ...")
out <- attach_to_obj(obj, out)

log$info("Saving results ...")
save_obj(out, outfile)
