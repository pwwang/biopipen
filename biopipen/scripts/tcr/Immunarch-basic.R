# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

log_info("")
log_info("# Basic analysis")
log_info("-----------------------------------")

# Fill up cases
fill_up_cases_basic = function(config) {
    log_debug("Filling up cases ...")
    cases = config$cases
    if (is.null(cases) || length(cases) == 0) {
        cases$DEFAULT = list(by = config$by, devpars = config$devpars, subset = config$subset)
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
            if (is.null(cases[[case]]$subset)) {
                cases[[case]]$subset = config$subset
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

do_one_case_basic = function(name, case, method) {
    log_info("- Processing case: {name} ...")
    odir = file.path(outdir, method)
    dir.create(odir, showWarnings = FALSE)

    if (!is.null(case$subset)) {
        d = immdata_from_expanded(filter_expanded_immdata(exdata, case$subset))
    } else {
        d = immdata
    }

    if (method == "len") {
        exp = repExplore(d$data, .method = method, .col = "aa")
    } else {
        exp = repExplore(d$data, .method = method)
    }
    if (is.null(case$by)) {
        p = vis(exp)
    } else {
        p = vis(exp, .by = case$by, .meta = d$meta)
    }

    ofig = file.path(odir, paste0(name, ".png"))
    png(ofig, width = case$devpars$width, height = case$devpars$height, res = case$devpars$res)
    print(p + scale_fill_biopipen())
    dev.off()

    ofig_pdf = file.path(odir, paste0(name, ".pdf"))
    pdf(ofig_pdf, width = case$devpars$width / case$devpars$res, height = case$devpars$height / case$devpars$res)
    print(p + scale_fill_biopipen())
    dev.off()

    add_report(
        list(
            src = ofig,
            name = if (name == "DEFAULT") NULL else name,
            download = ofig_pdf
        ),
        h1 = "Exploratory Analysis",
        h2 = switch(method,
            len = "CDR3 length distribution",
            volume = "Clonotype volume (Number of clonotypes)",
            count = "Clonotype abundances"
        ),
        ui = "table_of_images"
    )
}

# Do cases
do_cases_basic = function(cases, method) {
    log_info("Handling cases (method={method}) ...")
    for (name in names(cases)) {
        do_one_case_basic(name, cases[[name]], method)
    }
}

lens = fill_up_cases_basic(lens)
do_cases_basic(lens, "len")

volumes = fill_up_cases_basic(volumes)
do_cases_basic(volumes, "volume")

counts = fill_up_cases_basic(counts)
do_cases_basic(counts, "count")
