
# # https://stackoverflow.com/questions/50145643/unable-to-change-python-path-in-reticulate
# python = Sys.which({{envs.python | r}})
# Sys.setenv(RETICULATE_PYTHON = python)
# library(reticulate)
{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "single_cell.R" | source_r }}

library(immunarch)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)

immfile = {{in.immfile | r}}
outdir = normalizePath({{job.outdir | r}})
outfile = {{out.immfile | r}}
clusterfile = {{out.clusterfile | r}}
tool = {{envs.tool | r}}
python = {{envs.python | r}}
on_multi = {{envs.on_multi | r}}
args = {{envs.args | r}}
prefix = {{envs.prefix | r}}

setwd(outdir)

immdata = readRDS(immfile)
if (on_multi) {
    seqdata = immdata$multi
} else {
    seqdata = immdata$data
}
if (is.null(prefix)) { prefix = immdata$prefix }
if (is.null(prefix)) { prefix = "" }

get_cdr3aa_df = function() {
    out = expand_immdata(immdata, cell_id = "Barcode") %>%
        mutate(Barcode = glue(paste0(prefix, "{Barcode}")))

    if (on_multi) {
        out$CDR3.aa = sub(";", "", out$CDR3.aa)
    } else if ("chain" %in% colnames(out)) {
        out = out %>% separate_rows(chain, CDR3.aa, sep = ";") %>%
            filter(chain == "TRB")
    }
    out %>% select(Barcode, CDR3.aa)
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

clean_clustcr_output = function(clustcr_outfile, clustcr_input) {
    clustcr_out = read.delim2(clustcr_outfile, header=TRUE, row.names = NULL)
    colnames(clustcr_out) = c("CDR3.aa", "TCR_Cluster")
    in_cdr3 = read.delim2(clustcr_input, header=TRUE, row.names = NULL)
    out = left_join(in_cdr3, distinct(clustcr_out), by=c("CDR3.aa")) %>%
        mutate(
            TCR_Cluster = if_else(
                is.na(TCR_Cluster),
                paste0("S_", row_number()),
                paste0("M_", as.character(TCR_Cluster))
            )
        )
    out = left_join(
        cdr3aa_df,
        out,
        by = "CDR3.aa"
    )
    df = out %>%
        select(Barcode, TCR_Cluster) %>%
        add_count(TCR_Cluster, name="TCR_Cluster_Size") %>%
        distinct(Barcode, .keep_all = TRUE) %>%
        add_count(TCR_Cluster, name="TCR_Cluster_Size1") %>%
        column_to_rownames("Barcode")

    write.table(df, clusterfile, row.names=T, col.names=T, quote=F, sep="\t")
    out
}

run_clustcr = function() {
    log_info("Running ClusTCR ...")
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
    log_debug("- Running command: {clustcr_cmd}")
    rc = system(clustcr_cmd)
    if (rc != 0) {
        quit(status=rc)
    }
    clustcr_outfile = file.path(clustcr_dir, "clusters.txt")
    clean_clustcr_output(clustcr_outfile, clustcr_input)
}

prepare_giana = function() {
    giana_srcdir = "{{biopipen_dir}}/scripts/tcr/GIANA"

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
    # prepare input file for GIANA
    cdr3 = c()
    # cdr3col = if (!on_multi) "cdr3" else "CDR3.aa"
    cdr3col = "CDR3.aa"
    for (sample in names(seqdata)) {
        sdata = seqdata[[sample]]
        if (on_multi) {
            sdata[[cdr3col]] = sub(";", "", sdata[[cdr3col]])
        } else if ("chain" %in% colnames(sdata)) {
            sdata = sdata %>% separate_rows(chain, cdr3col, sep = ";") %>%
                filter(chain == "TRB")
        }
        cdr3 = union(cdr3, unique(sdata[[cdr3col]]))
    }
    cdr3 = unique(cdr3)

    # cdr3 = distinct(cdr3, aminoAcid, vMaxResolved)

    cdr3file = file.path(outdir, "cdr3.csv")
    write.table(
        data.frame(CDR3.aa=cdr3),
        cdr3file,
        row.names=FALSE, col.names=TRUE, quote=FALSE
    )
    cdr3file
}

clean_giana_output = function(giana_outfile, giana_infile) {
    # generate an output file with columns:
    # CDR3.aa, TCR_Cluster, V.name, Sample
    # If sequence doesn't exist in the input file,
    # Then a unique cluster id is assigned to it.
    giana_out = read.delim2(giana_outfile, header=FALSE, comment.char = "#", row.names = NULL)[, 1:2, drop=FALSE]
    colnames(giana_out) = c("CDR3.aa", "TCR_Cluster")
    in_cdr3 = read.delim2(giana_infile, header=TRUE, row.names = NULL)
    out = left_join(in_cdr3, distinct(giana_out), by=c("CDR3.aa")) %>%
        mutate(
            TCR_Cluster = if_else(
                is.na(TCR_Cluster),
                paste0("S_", row_number()),
                paste0("M_", as.character(TCR_Cluster))
            )
        )

    out = left_join(
        cdr3aa_df,
        out,
        by = "CDR3.aa"
    )
    df = out %>%
        select(Barcode, TCR_Cluster) %>%
        add_count(TCR_Cluster, name="TCR_Cluster_Size") %>%
        distinct(Barcode, .keep_all = TRUE) %>%
        add_count(TCR_Cluster, name="TCR_Cluster_Size1") %>%
        column_to_rownames("Barcode")

    write.table(df, clusterfile, row.names=T, col.names=T, quote=F, sep="\t")
    out
}

run_giana = function() {
    log_info("Running GIANA ...")
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
    log_debug("- Running command: {giana_cmd}")
    rc = system(giana_cmd)
    if (rc != 0) {
        quit(status=rc)
    }
    giana_outfile = file.path(giana_outdir, "cdr3--RotationEncodingBL62.txt")
    clean_giana_output(giana_outfile, giana_input)
}

attach_to_immdata = function(out) {
    seqdata2 = list()
    # by = if (!on_multi) c(cdr3 = "CDR3.aa") else "CDR3.aa"
    by = "CDR3.aa"
    for (sample in names(seqdata)) {
        sample_out = left_join(seqdata[[sample]], out, by=by)
        seqdata2[[sample]] = sample_out
        if (!on_multi) {
            immdata$data[[sample]] = immdata$data[[sample]] %>% left_join(
                out, by = "CDR3.aa"
            )
        } else {
            immdata$multi[[sample]] = immdata$multi[[sample]] %>% left_join(
                out, by = c(cdr3 = "CDR3.aa")
            )
        }
        # if ("single" %in% names(immdata)) {
        #     immdata$data[[sample]] = immdata$data[[sample]] %>% left_join(
        #         out, by = "CDR3.aa"
        #     )
        # }
    }
    if (!on_multi) {
        immdata$data = seqdata2
    } else {
        immdata$multi = seqdata2
    }
    saveRDS(immdata, file = outfile)
    # seqdata2
}


if (tolower(tool) == "clustcr") {
    out = run_clustcr()
} else if (tolower(tool) == "giana") {
    out = run_giana()
} else {
    stop(paste("Unknown tool:", tool))
}

log_info("Saving results ...")
attach_to_immdata(out)
