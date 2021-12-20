library(Seurat)
library(dplyr)

srtobjfile = {{in.srtobj | r}}
out = {{out.out | r}}
invert = {{envs.invert | r}}
filterfmt = {{envs.filterfmt | r}}
multicase = {{envs.multicase | r}}

srtobj = readRDS(srtobjfile)

{% set filtercontent = in.filterfile | read %}
{% set filterfmt = envs.filterfmt %}
{% if envs.filterfmt == "auto" %}
{%    if "=" in filtercontent %}
{%         set filterfmt = "subset" %}
{%    else %}
{%         set filterfmt = "grouping" %}
{%    endif %}
{% endif %}

{% if filterfmt == "subset" %}

    {% set subsets = in.filterfile | read | toml_loads %}
    {% for key, val in subsets.subsetting.items() %}
        sobj = subset(srtobj, subset = {{val}}, invert = invert)
        if (multicase) {
            outfile = file.path(out, "{{key}}.RDS")
        } else {
            outfile = out
        }
        saveRDS(sobj, outfile)
    {% endfor %}

{% else %}

    groups = read.table(
        {{ in.filterfile | r }},
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
    } else {
        groups = groups %>%
            mutate(ALL=strsplit(ALL, ";", fixed=TRUE))
    }

    subsetSeurat_group_single = function(srtobj, name, cells) {
        sobj = subset(srtobj, cells = cells, invert = invert)
        if (multicase) {
            outfile = file.path(out, paste0(name, ".RDS"))
        } else {
            outfile = out
        }
        saveRDS(sobj, outfile)
    }

    subsetSeurat_group = function(srtobj, groups) {
        for (i in seq_len(nrow(groups))) {
            row = groups[i,,drop=TRUE]
            name = as.character(row[1])
            cells = unlist(row$ALL)
            subsetSeurat_group_single(srtobj, name, cells)
        }
    }

    subsetSeurat_group(srtobj, groups)

{% endif %}
