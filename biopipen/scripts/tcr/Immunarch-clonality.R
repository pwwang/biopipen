# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

top_clones = {{envs.top_clones | r}}
rare_clones = {{envs.rare_clones | r}}
hom_clones = {{envs.hom_clones | r}}

# Fill up cases
fill_up_cases = function(config) {
    cases = config$cases
    if (is.null(cases) || length(cases) == 0) {
        cases$DEFAULT = list(by = config$by, marks = config$marks, devpars = config$devpars)
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
    odir = file.path(outdir, paste0(method, "_clones"))
    dir.create(odir, showWarnings = FALSE)
    if (method == "top") {
        exp = repClonality(immdata$data, .method = method, .head = case$marks)
    } else if (method == "rare") {
        exp = repClonality(immdata$data, .method = method, .bound = case$marks)
    } else if (method == "homeo") {
        marks = case$marks
        if (is.null(marks$Rare)) { marks$Rare = 1e-5 }
        if (is.null(marks$Small)) { marks$Small = 1e-4 }
        if (is.null(marks$Medium)) { marks$Medium = 1e-3 }
        if (is.null(marks$Large)) { marks$Large = 1e-2 }
        if (is.null(marks$Hyperexpanded)) { marks$Hyperexpanded = 1 }
        marks = marks[c("Rare", "Small", "Medium", "Large", "Hyperexpanded")]
        exp = repClonality(immdata$data, .method = method, .clone.types = unlist(marks))
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
    print(paste0("- Clonality analysls: ", method))
    for (name in names(cases)) {
        do_one_case(name, cases[[name]], method)
    }
}

top_clones = fill_up_cases(top_clones)
do_cases(top_clones, "top")

rare_clones = fill_up_cases(rare_clones)
do_cases(rare_clones, "rare")

hom_clones = fill_up_cases(hom_clones)
do_cases(hom_clones, "homeo")
