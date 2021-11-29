library(dplyr)
library(immunarch)

immfile = {{in.immdata | quote}}
config = {{in.configfile | read | toml_loads | r}}
outfile = {{out.outfile | quote}}
groupfile = {{out.groupfile | quote}}

immdata0 = readRDS(immfile)

groupdata = NULL
for (name in names(config)) {
    conf = config[[name]]
    clonotype = if (is.null(conf$clonotype)) FALSE else conf$clonotype
    filterdata = if (is.null(conf$filterdata)) TRUE else conf$filterdata
    filters = if (is.null(conf$filters)) list() else conf$filters
    immdata = immdata0

    filters_order = lapply(
        filters,
        function(x) if(is.null(x$ORDER)) 0 else x$ORDER
    )
    by_order = names(filters_order[order(unlist(filters_order))])

    for (by in by_order) {
        cond = filters[[by]]
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

    if (filterdata) {
        saveRDS(immdata, file=outfile)
    }

    # Save a group file
    if (!clonotype) {
        groupdf = list()
        for (sample in names(immdata$data)) {
            groupdf[[sample]] = paste(
                immdata$data[[sample]]$Barcode,
                collapse=";"
            )
        }
        cnames = names(groupdf)
        groupdf = as.data.frame(groupdf)
        colnames(groupdf) = cnames
        groupdf$Group = name
        groupdf = groupdf %>% select(Group, everything())
    } else {
        groupdf = NULL
        for (sample in names(immdata$data)) {
            df = immdata$data[[sample]] %>% select(Clonotype=CDR3.aa, Barcode)
            colnames(df)[2] = sample
            if (is.null(groupdf)) {
                groupdf = df
            } else {
                groupdf = full_join(groupdf, df, by="Clonotype")
            }
        }
    }
    if (is.null(groupdata)) {
        groupdata = groupdf
    } else if (colnames(groupdata) != colnames(groupdf)) {
        stop("Samples are not the same after filtering")
    } else {
        groupdata = bind_rows(groupdata, groupdf)
    }
}

write.table(groupdf, groupfile, quote=F, row.names=F, col.names=T, sep="\t")
