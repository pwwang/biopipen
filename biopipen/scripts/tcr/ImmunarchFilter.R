library(dplyr)
library(immunarch)

immfile = {{in.immdata | quote}}
outfile = {{out.outfile | quote}}
groupfile = {{out.groupfile | quote}}
condition = {{in.filters | read | toml_loads | r}}
groupname = {{(in.name if in.name else "Filtered") | quote}}

if (endsWith(groupname, "@clonotype")) {
    groupname = substr(groupname, 1, nchar(groupname) - 10)
    by_clonotype = TRUE
} else {
    by_clonotype = FALSE
}

condition_order = lapply(condition, function(x) ifelse(is.null(x$ORDER), 0, x$ORDER))
by_order = names(condition_order[order(unlist(condition_order))])

if (length(condition) == 0) {
    stop("No condition provided.")
}

immdata = readRDS(immfile)

for (by in by_order) {
    cond = condition[[by]]
    cond$ORDER = NULL
    if (by == "by.count") {

        cdr3aa = c()
        for (sample in names(immdata$data)) {
            cdr3aa = unique(c(cdr3aa, immdata$data[[sample]]$CDR3.aa))
        }
        counts = data.frame(CDR3.aa = cdr3aa)
        for (sample in names(immdata$data)) {
            counts = left_join(counts, immdata$data[[sample]], by="CDR3.aa") %>%
                select(all_of(colnames(counts)), Clones)
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

saveRDS(immdata, file=outfile)

# Save a group file
if (!by_clonotype) {
    groupdata = list()
    for (sample in names(immdata$data)) {
        groupdata[[sample]] = paste(immdata$data[[sample]]$Barcode, collapse=";")
    }
    groupdf = as.data.frame(groupdata)
    colnames(groupdf) = names(groupdata)
    groupdf$Group = groupname
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

write.table(groupdf, groupfile, quote=F, row.names=F, col.names=T, sep="\t")
