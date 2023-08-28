# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

volumes = {{envs.volumes | r}}
lens = {{envs.lens | r}}
counts = {{envs.counts | r}}

# Fill up cases
fill_up_cases = function(config) {
    cases = config$cases
    if (is.null(cases) || length(cases) == 0) {
        cases$DEFAULT = list(by = config$by, devpars = config$devpars)
    } else {
        for (case in names(cases)) {
            if (is.null(cases[[case]]$by)) {
                cases[[case]]$by = config$by
            }
            if (is.null(cases[[case]]$devpars)) {
                cases[[case]]$devpars = config$devpars
            }
            if (is.null(cases[[case]]$devpars$width)) {
                cases[[case]]$devpars$width = config$devpars$width
            }
            if (is.null(cases[[case]]$devpars$height)) {
                cases[[case]]$devpars$height = config$devpars$height
            }
            if (is.null(cases[[case]]$devpars$res)) {
                cases[[case]]$devpars$res = config$devpars$res
            }
        }
    }
    for (case in names(cases)) {
        by = cases[[case]]$by
        by = if (is.character(by)) trimws(unlist(strsplit(by, ","))) else NULL
        cases[[case]]$by = by
    }

    cases
}

do_one_case = function(name, case, method) {
    print(paste0("  Case: ", name))
    odir = file.path(outdir, method)
    dir.create(odir, showWarnings = FALSE)
    if (method == "len") {
        exp = repExplore(immdata$data, .method = method, .col = "aa")
    } else {
        exp = repExplore(immdata$data, .method = method)
    }
    if (is.null(case$by)) {
        p = vis(exp)
    } else {
        p = vis(exp, .by = case$by, .meta = immdata$meta)
    }
    ofig = file.path(odir, paste0(name, ".png"))
    png(ofig, width = case$devpars$width, height = case$devpars$height, res = case$devpars$res)
    print(p)
    dev.off()
}

# Do cases
do_cases = function(cases, method) {
    print(paste0("- Basic analysls: ", method))
    for (name in names(cases)) {
        do_one_case(name, cases[[name]], method)
    }
}

volumes = fill_up_cases(volumes)
do_cases(volumes, "volume")

lens = fill_up_cases(lens)
do_cases(lens, "len")

counts = fill_up_cases(counts)
do_cases(counts, "count")
