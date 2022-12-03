library(Seurat)
library(rlang)
library(dplyr)

srtobjfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
ignore_nas = {{envs.ignore_nas | r}}

srtobj = readRDS(srtobjfile)

{% set configs = in.subsets | config: "toml" %}
{% for casename, caseconf in configs.items() %}
    casename = {{ casename | r }}
    mutaters = {{ caseconf | attr: "mutaters" | r }}
    groupby = {{ caseconf | attr: "groupby" | r }}
    subset_cond = {{ caseconf | attr: "subset" | r }}
    print(paste("- Working on case:", casename))
    srtobj_copy = srtobj
    if (!is.null(mutaters)) {
        expr = list()
        for (key in names(mutaters)) {
            expr[[key]] = parse_expr(mutaters[[key]])
        }
        srtobj_copy@meta.data = srtobj_copy@meta.data %>% mutate(!!!expr)
    }

    if (!is.null(groupby)) {
        print(paste("  with groupby:", groupby))
        if (ignore_nas) {
            cells = FetchData(srtobj_copy, vars = groupby) %>%
                filter(!is.na(.data[[groupby]])) %>%
                rownames()

            srtobj_copy = subset(srtobj_copy, cells = cells)
        } else {
            # convert NA to string
            srtobj_copy@meta.data[[groupby]] = as.character(srtobj_copy@meta.data[[groupby]])
        }
        sobjs = SplitObject(srtobj_copy, split.by = groupby)
        for (cname in names(sobjs)) {
            outfile = file.path(outdir, paste0(casename, "_", cname, ".RDS"))
            saveRDS(sobjs[[cname]], outfile)
        }
    } else {
        print(paste("  with subsetting:", subset_cond))
        outfile = file.path(outdir, paste0(casename, ".RDS"))
        subset_cmd = paste0("subset(srtobj_copy, subset = ", subset_cond, ")")
        sobj = eval(parse(text = subset_cmd))
        saveRDS(sobj, outfile)
    }

{% endfor %}
