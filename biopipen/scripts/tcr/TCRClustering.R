
# # https://stackoverflow.com/questions/50145643/unable-to-change-python-path-in-reticulate
# python = Sys.which({{envs.python | r}})
# Sys.setenv(RETICULATE_PYTHON = python)
# library(reticulate)

library(ComplexHeatmap)
library(immunarch)
library(dplyr)
library(tidyr)

immfile = {{in.immfile | r}}
outdir = {{job.outdir | r}}
outfile = {{out.immfile | r}}
clusterfile = {{out.clusterfile | r}}
heatmap = {{out.heatmap | r}}
tool = {{envs.tool | r}}
python = {{envs.python | r}}
tmpdir = {{envs.tmpdir | r}}
giana_source = {{envs.giana_source | r}}
args = {{envs.args | r}}
heatmap_meta = {{envs.heatmap_meta | r}}

immdata = readRDS(immfile)
if ("single" %in% names(immdata)) {
    seqdata = immdata$single
} else {
    seqdata = immdata$data
}

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
cdr3df = pd.read_csv(clustcr_infile, sep="\t", index_col=None)
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
    colnames(clustcr_out) = c(
        "CDR3.aa",
        "TCR_Cluster"
    )
    in_cdr3 = read.delim2(clustcr_input, header=TRUE, row.names = NULL)
    colnames(in_cdr3) = c("CDR3.aa", "V.name", "Sample")
    out = left_join(in_cdr3, distinct(clustcr_out), by=c("CDR3.aa")) %>%
        mutate(
            V.name = if_else(
                endsWith(V.name, "*01"),
                substr(V.name, 1, nchar(V.name) - 3),
                V.name
            ),
            TCR_Cluster = if_else(
                is.na(TCR_Cluster),
                paste0("UNK_", row_number()),
                as.character(TCR_Cluster)
            )
        )
    write.table(out, clusterfile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    out
}

run_clustcr = function() {
    print(paste("Using tool:", "ClusTCR"))
    clustcr_dir = file.path(outdir, "ClusTCR_Output")
    dir.create(clustcr_dir, showWarnings = FALSE)
    clustcr_file = prepare_clustcr(clustcr_dir)
    clustcr_input = prepare_giana_input()
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
        download.file(giana_source$url, giana_file)
        download.file(giana_source$giana4, giana4_file)
        download.file(giana_source$query, giana_query)
        download.file(giana_source$trbv, giana_trbv)
    }

    giana_srcdir
}

prepare_giana_input = function() {
    # prepare input file for GIANA
    cdr3 = c()
    for (sample in names(seqdata)) {
        cdr3 = bind_rows(cdr3, seqdata[[sample]] %>%
            transmute(aminoAcid=CDR3.aa, vMaxResolved=paste0(V.name, "*01"), Sample=sample))
    }

    # cdr3 = distinct(cdr3, aminoAcid, vMaxResolved)

    cdr3file = file.path(outdir, "cdr3.csv")
    write.table(cdr3, cdr3file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    cdr3file
}

clean_giana_output = function(giana_outfile, giana_infile) {
    # generate an output file with columns:
    # CDR3.aa, TCR_Cluster, V.name, Sample
    # If sequence doesn't exist in the input file,
    # Then a unique cluster id is assigned to it.
    giana_out = read.delim2(giana_outfile, header=FALSE, comment.char = "#", row.names = NULL)
    colnames(giana_out) = c(
        "CDR3.aa",
        "TCR_Cluster",
        "V.name",
        "Sample"
    )
    in_cdr3 = read.delim2(giana_infile, header=TRUE, row.names = NULL)
    colnames(in_cdr3) = c("CDR3.aa", "V.name", "Sample")
    out = left_join(in_cdr3, distinct(giana_out), by=c("CDR3.aa", "V.name", "Sample")) %>%
        mutate(
            V.name = if_else(
                endsWith(V.name, "*01"),
                substr(V.name, 1, nchar(V.name) - 3),
                V.name
            ),
            TCR_Cluster = if_else(
                is.na(TCR_Cluster),
                paste0("UNK_", row_number()),
                as.character(TCR_Cluster)
            )
        )
    write.table(out, clusterfile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    out
}

run_giana = function() {
    print(paste("Using tool:", "GIANA"))
    giana_srcdir = prepare_giana()
    giana_input = prepare_giana_input()
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
    for (sample in names(seqdata)) {
        sample_out = filter(out, Sample == sample) %>% select(-"Sample")
        sample_out = left_join(seqdata[[sample]], sample_out, by=c("CDR3.aa", "V.name"))
        seqdata2[[sample]] = sample_out
    }
    if ("single" %in% names(immdata)) {
        immdata$single = seqdata2
    } else {
        immdata$data = seqdata2
    }
    saveRDS(immdata, file = outfile)
}

plot_heatmap = function(out) {
    # Generate data for heatmap
    samples = names(seqdata)
    tcr_clusters = list()
    for (sample in samples) {
        tcr_clusters[[sample]] = filter(out, Sample == sample & !is.na(TCR_Cluster)) %>%
            pull(TCR_Cluster) %>% unique()
    }
    plotdata = matrix(NA, ncol = length(samples), nrow = length(samples))
    rownames(plotdata) = samples
    colnames(plotdata) = samples
    for (sample1 in samples) {
        for (sample2 in samples) {
            plotdata[sample1, sample2] = length(
                intersect(tcr_clusters[[sample1]], tcr_clusters[[sample2]])
            )
        }
    }

    if (length(heatmap_meta) == 0) {
        anno = NULL
    } else {
        anno = as.list(immdata$meta[, heatmap_meta, drop=FALSE])
        anno = do.call(HeatmapAnnotation, anno)
    }

    # Plot heatmap
    hm = Heatmap(
        plotdata,
        name = "Shared TCR Clusters",
        col = c("#ffe1e1", "red3"),
        cluster_rows = FALSE,
        top_annotation = anno,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(plotdata[samples[i], samples[j]], x, y, gp = gpar(fontsize = 10))
        }
    )
    png(heatmap, width = if (is.null(anno)) 1000 else 1200, height = 1000, res = 100)
    print(hm)
    dev.off()

}

if (tolower(tool) == "clustcr") {
    out = run_clustcr()
} else if (tolower(tool) == "giana") {
    out = run_giana()
} else {
    stop(paste("Unknown tool:", tool))
}

attach_to_immdata(out)
plot_heatmap(out)
