# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

log_info("")
log_info("# Clonality analysis")
log_info("-----------------------------------")

# Fill up cases
fill_up_cases_clonality = function(config) {
    cases = config$cases
    if (is.null(cases) || length(cases) == 0) {
        cases$DEFAULT = list(by = config$by, marks = config$marks, devpars = config$devpars, subset = config$subset)
    } else {
        for (case in names(cases)) {
            if (is.null(cases[[case]]$by)) {
                cases[[case]]$by = config$by
            }
            if (is.null(cases[[case]]$marks)) {
                cases[[case]]$marks = config$marks
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

do_one_case_clonality = function(name, case, method) {
    log_info("- Processing case: {name} ...")
    odir = file.path(outdir, paste0(method, "_clones"))
    dir.create(odir, showWarnings = FALSE)

    if (is.null(case$subset)) {
        d = immdata
    } else {
        d = immdata_from_expanded(filter_expanded_immdata(exdata, case$subset))
    }
    if (method == "top") {
        exp = repClonality(d$data, .method = method, .head = case$marks)
    } else if (method == "rare") {
        exp = repClonality(d$data, .method = method, .bound = case$marks)
    } else if (method == "homeo") {
        marks = case$marks
        if (is.null(marks$Rare)) { marks$Rare = 1e-5 }
        if (is.null(marks$Small)) { marks$Small = 1e-4 }
        if (is.null(marks$Medium)) { marks$Medium = 1e-3 }
        if (is.null(marks$Large)) { marks$Large = 1e-2 }
        if (is.null(marks$Hyperexpanded)) { marks$Hyperexpanded = 1 }
        marks = marks[c("Rare", "Small", "Medium", "Large", "Hyperexpanded")]
        exp = repClonality(d$data, .method = method, .clone.types = unlist(marks))
    }
    if (is.null(case$by)) {
        p = vis(exp)
    } else {
        p = vis(exp, .by = case$by, .meta = d$meta)
    }
    ofig = file.path(odir, paste0(name, ".png"))
    png(ofig, width = case$devpars$width, height = case$devpars$height, res = case$devpars$res)
    print(p)
    dev.off()

    ofig_pdf = file.path(odir, paste0(name, ".pdf"))
    pdf(ofig_pdf, width = case$devpars$width / case$devpars$res, height = case$devpars$height / case$devpars$res)
    print(p)
    dev.off()

    add_report(
        list(
            src = ofig,
            name = if (name == "DEFAULT") NULL else name
        ),
        h1 = "Clonality Analysis",
        h2 = switch(method,
            top = "Top Clones",
            rare = "Rare Clones",
            homeo = "Clonal Space Homeostasis"
        ),
        ui = "table_of_images"
    )
}

# Do cases
do_cases_clonality = function(cases, method) {
    log_info("Handling cases (method={method}) ...")
    for (name in names(cases)) {
        do_one_case_clonality(name, cases[[name]], method)
    }
}

top_clones = fill_up_cases_clonality(top_clones)
do_cases_clonality(top_clones, "top")

rare_clones = fill_up_cases_clonality(rare_clones)
do_cases_clonality(rare_clones, "rare")

hom_clones = fill_up_cases_clonality(hom_clones)
do_cases_clonality(hom_clones, "homeo")
