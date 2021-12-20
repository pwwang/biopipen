library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)

srtdir = {{in.srtdir | r}}
groupfile = {{in.groupfile | r}}
config = {{in.configfile | read | toml_loads | r}}
outdir = {{out.outdir | r}}

do_one_srtobjfile = function(srtobjfile) {
    srtobj = readRDS(srtobjfile)
    outfile = file.path(outdir, basename(srtobjfile))

    groupCellsByInput = function(srtobj, groupfile) {
        groups = read.table(
            groupfile,
            sep = "\t",
            header = TRUE,
            row.names = NULL,
            check.names = FALSE
        )
        cnames = colnames(groups)
        gname = cnames[1]
        samples = cnames[2:length(cnames)]
        if (samples != "ALL") {
            groups = groups %>%
                group_by(across(all_of(gname))) %>%
                mutate(across(all_of(samples), ~ strsplit(.x, ";", fixed=TRUE))) %>%
                summarise(across(all_of(samples), ~ list(c(unlist(.x))))) %>%
                rowwise() %>%
                mutate(
                    across(
                        all_of(samples),
                        ~ list(paste(cur_column(), unlist(.x), sep="_"))
                    )
                ) %>%
                mutate(ALL=list(c_across(all_of(samples)))) %>%
                select(all_of(gname), "ALL")
        }

        metadata = srtobj$meta.data
        metadata[[config$grouping_name]] = ""
        for (i in seq_len(nrow(groups))) {
            row = groups[i,,drop=TRUE]
            metadata[unlist(row$ALL), config$grouping_name] = row[[gname]]
        }
        srtobj$meta.data = metadata

        saveRDS(srtobj, file = outfile)
    }


    groupCellsByIdents = function(srtobj) {
        srtobj[[config$grouping_name]] = Idents(srtobj)
        saveRDS(srtobj, file = outfile)
    }


    groupCellsByConfig = function(srtobj, groups) {
        metadata = srtobj$meta.data
        case_when_conds = with(metadata, {
            lapply(names(groups), function(name) {
                eval(parse(text = groups[[name]])) ~ name
            })
        })
        cells = rownames(metadata)
        metadata = metadata %>%
            mutate(do.call(case_when, case_when_cond)) %>%
            as.data.frame()
        colnames(metadata)[ncol(metadata)] = config$grouping_name
        rownames(metadata) = cells
        srtobj$meta.data = metadata

        saveRDS(srtobj, file = outfile)
    }


    {% set grouping = in.configfile | read | toml_loads | attr: "grouping" %}
    {% if grouping == "Input" %}
    groupCellsByInput(srtobj, groupfile)
    {% elif grouping == "Idents" %}
    groupCellsByIdents(srtobj)
    {% else %}
    groupCellsByConfig(
        srtobj,
        {{in.groupfile | read | toml_loads | r}}
    )
    {% endif %}
}

for (srtobjfile in Sys.glob(file.path(srtdir, "*.RDS"))) {
    do_one_srtobjfile(srtobjfile)
}
