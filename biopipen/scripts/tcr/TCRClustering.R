
# # https://stackoverflow.com/questions/50145643/unable-to-change-python-path-in-reticulate
# python = Sys.which({{envs.python | r}})
# Sys.setenv(RETICULATE_PYTHON = python)
# library(reticulate)

library(immunarch)
library(dplyr)
library(tidyr)
library(tibble)

immfile = {{in.immfile | r}}
outdir = normalizePath({{job.outdir | r}})
outfile = {{out.immfile | r}}
clusterfile = {{out.clusterfile | r}}
tool = {{envs.tool | r}}
python = {{envs.python | r}}
tmpdir = {{envs.tmpdir | r}}
on_multi = {{envs.on_multi | r}}
giana_repo = {{envs.giana_repo | r}}
args = {{envs.args | r}}

setwd(outdir)

immdata = readRDS(immfile)
if (on_multi) {
    seqdata = immdata$multi
} else {
    seqdata = immdata$data
}

get_cdr3aa_df = function() {
    out = NULL
    for (sample in names(immdata$data)) {
        tmpdf = immdata$data[[sample]] |>
            select(Barcode, CDR3.aa) |>
            separate_rows(Barcode, sep = ";") |>
            mutate(Barcode = paste0(sample, "_", Barcode))
        out = bind_rows(out, tmpdf)
    }
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
import pandas as pd
import clustcr

clustcr_dir, clustcr_infile = sys.argv[1:3]
cdr3df = pd.read_csv(clustcr_infile, index_col=None)
cdr3 = cdr3df.iloc[:, 0]

clustering = clustcr.Clustering(%s)
output = clustering.fit(cdr3)
output.clusters_df.to_csv(clustcr_dir + "/clusters.txt", sep="\\t", index=False)
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
    df = out |>
        select(Barcode, TCR_Cluster) |>
        distinct(Barcode, .keep_all = TRUE) |>
        column_to_rownames("Barcode")

    write.table(df, clusterfile, row.names=T, col.names=T, quote=F, sep="\t")
    out
}

run_clustcr = function() {
    print(paste("Using tool:", "ClusTCR"))
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
    rc = system(clustcr_cmd)
    if (rc != 0) {
        quit(status=rc)
    }
    clustcr_outfile = file.path(clustcr_dir, "clusters.txt")
    clean_clustcr_output(clustcr_outfile, clustcr_input)
}

prepare_giana = function() {
    giana_srcdir = file.path(tmpdir, "GIANA_SOURCE")
    dir.create(giana_srcdir, showWarnings = FALSE)

    giana_file = file.path(giana_srcdir, "GIANA.py")
    giana4_file = file.path(giana_srcdir, "GIANA4.py")
    giana_query = file.path(giana_srcdir, "query.py")
    giana_trbv = file.path(giana_srcdir, "Imgt_Human_TRBV.fasta")
    if (!file.exists(giana_file)) {
        download.file(paste(giana_repo, "GIANA4.1.py", sep="/"), giana_file)
        download.file(paste(giana_repo, "GIANA4.py", sep="/"), giana4_file)
        download.file(paste(giana_repo, "query.py", sep="/"), giana_query)
        download.file(paste(giana_repo, "Imgt_Human_TRBV.fasta", sep="/"), giana_trbv)
    }

    giana_srcdir
}

prepare_input = function() {
    # prepare input file for GIANA
    cdr3 = c()
    # cdr3col = if (!on_multi) "cdr3" else "CDR3.aa"
    cdr3col = "CDR3.aa"
    for (sample in names(seqdata)) {
        # cdr3 = bind_rows(cdr3, seqdata[[sample]] %>%
        #     transmute(aminoAcid=CDR3.aa, vMaxResolved=paste0(V.name, "*01"), Sample=sample))
        cdr3 = union(
            cdr3,
            seqdata[[sample]] %>% pull(cdr3col) %>% unique()
        )
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
    df = out |>
        select(Barcode, TCR_Cluster) |>
        distinct(Barcode, .keep_all = TRUE) |>
        column_to_rownames("Barcode")

    write.table(df, clusterfile, row.names=T, col.names=T, quote=F, sep="\t")
    out
}

run_giana = function() {
    print(paste("Using tool:", "GIANA"))
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

attach_to_immdata(out)
