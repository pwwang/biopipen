{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(dplyr)
library(glue)
library(tidyr)
library(tibble)
library(immunarch)

immfile = {{in.immdata | r}}
{% if in.filterfile %}
filters = {{in.filterfile | toml_load | r}}
{% else %}
filters = {{envs.filters | r}}
{% endif %}
metacols = {{envs.metacols | r}}
outfile = {{out.outfile | r}}
groupfile = {{out.groupfile | r}}

immdata0 = readRDS(immfile)
groupname = filters$name

groupdata = NULL
for (name in names(filters$filters)) {
    filt = filters$filters[[name]]
    if (is.null(filt)) {
        rest = name
        next
    }

    immdata = immdata0
    filters_order = lapply(
        filt,
        function(x) if(is.null(x$ORDER)) 0 else x$ORDER
    )
    by_order = names(filters_order[order(unlist(filters_order))])
    save = if (is.null(filt$SAVE)) FALSE else filt$SAVE
    filt$SAVE = NULL

    for (by in by_order) {
        cond = filt[[by]]
        cond$ORDER = NULL
        if (by == "by.count") {
            cdr3aa = c()
            for (sample in names(immdata$data)) {
                cdr3aa = unique(c(cdr3aa, immdata$data[[sample]]$CDR3.aa))
            }
            counts = data.frame(CDR3.aa = cdr3aa)
            for (sample in names(immdata$data)) {
                counts = left_join(
                    counts,
                    immdata$data[[sample]],
                    by="CDR3.aa"
                ) %>% select(all_of(colnames(counts)), Clones)
                colnames(counts)[ncol(counts)] = sample
            }
            counts = counts %>% distinct(CDR3.aa, .keep_all = TRUE)
            counts[is.na(counts)] = 0
            counts = counts %>%
                mutate(TOTAL = rowSums(select(., -"CDR3.aa"))) %>%
                arrange(desc(TOTAL)) %>%
                filter(eval(parse(text=cond$filter)))

            for (sample in names(immdata$data)) {
                immdata$data[[sample]] = immdata$data[[sample]] %>%
                    filter(CDR3.aa %in% counts$CDR3.aa)
            }
        } else {
            cond2 = list()
            for (condname in names(cond)) {
                con = cond[[condname]]
                if (startsWith(con, "include(") ||
                    startsWith(con, "exclude(") ||
                    startsWith(con, "morethan(") ||
                    startsWith(con, "lessthan(") ||
                    startsWith(con, "interval(")
                ) {
                    cond2[[condname]] = eval(parse(text=con))
                } else {
                    cond2[[condname]] = include(con)
                }
            }
            immdata = repFilter(
                immdata,
                .method = by,
                .query = cond2
            )
        }
    }

    if (save) {
        if (file.exists(outfile)) {
            stop("Can only save data with one filter.")
        }
        saveRDS(immdata, file=outfile)
    }

    # Save a group file
    groupdf = do_call(rbind, lapply(
        names(immdata$data),
        function(sample) {
            mdata = as.list(immdata$meta[immdata$meta$Sample == sample, ])
            data = immdata$data[[sample]] %>%
                separate_rows(Barcode, sep=";") %>%
                distinct(Barcode, .keep_all = TRUE)
            for (mname in names(mdata)) {
                data[[mname]] = mdata[[mname]]
                assign(mname, mdata[[mname]])
            }
            Barcode = data$Barcode

            gdf = data %>%
                select(all_of(metacols)) %>%
                mutate(Cell = glue("{{envs.prefix}}{Barcode}"))
            gdf[[groupname]] = name
            gdf
        }
    ))

    if (is.null(groupdata)) {
        groupdata = groupdf
    } else {
        groupdata = rbind(groupdata, groupdf)
    }
}

if (!file.exists(outfile)) {
    saveRDS(immdata0, file=outfile)
}
groupdata = groupdata %>%
    as.data.frame() %>%
    column_to_rownames("Cell")
write.table(groupdata, groupfile, quote=F, row.names=T, col.names=T, sep="\t")
