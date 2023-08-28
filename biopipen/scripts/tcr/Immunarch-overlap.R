# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

overlaps = {{ envs.overlaps | r: todot="-" }}

# Fill up cases
cases = overlaps$cases
if (is.null(cases) || length(cases) == 0) {
    cases[[overlaps$method]] = list(
        method = overlaps$method,
        vis_pars = overlaps$vis_pars,
        devpars = overlaps$devpars,
        analyses = overlaps$analyses
    )
} else {
    for (name in names(cases)) {
        if (is.null(cases[[name]]$method)) {
            cases[[name]]$method = overlaps$method
        }
        if (is.null(cases[[name]]$vis_pars)) {
            cases[[name]]$vis_pars = overlaps$vis_pars
        }
        if (is.null(cases[[name]]$devpars)) {
            cases[[name]]$devpars = overlaps$devpars
        }
        if (is.null(cases[[name]]$analyses)) {
            cases[[name]]$analyses = overlaps$analyses
        }
        if (is.null(cases[[name]]$devpars$width)) {
            cases[[name]]$devpars$width = overlaps$devpars$width
        }
        if (is.null(cases[[name]]$devpars$height)) {
            cases[[name]]$devpars$height = overlaps$devpars$height
        }
        if (is.null(cases[[name]]$devpars$res)) {
            cases[[name]]$devpars$res = overlaps$devpars$res
        }
    }
}

# analyses cases
for (name in names(cases)) {
    analyses = cases[[name]]$analyses
    if (is.null(analyses$cases) || length(analyses$cases) == 0) {
        analyses$cases$DEFAULT = list(
            method = analyses$method,
            vis_args = analyses$vis_args,
            devpars = analyses$devpars
        )
    } else {
        for (aname in names(analyses$cases)) {
            if (is.null(analyses$cases[[aname]]$method)) {
                analyses$cases[[aname]]$method = analyses$method
            }
            if (is.null(analyses$cases[[aname]]$vis_args)) {
                analyses$cases[[aname]]$vis_args = analyses$vis_args
            }
            if (is.null(analyses$cases[[aname]]$devpars)) {
                analyses$cases[[aname]]$devpars = analyses$devpars
            }
            if (is.null(analyses$cases[[aname]]$devpars$width)) {
                analyses$cases[[aname]]$devpars$width = analyses$devpars$width
            }
            if (is.null(analyses$cases[[aname]]$devpars$height)) {
                analyses$cases[[aname]]$devpars$height = analyses$devpars$height
            }
            if (is.null(analyses$cases[[aname]]$devpars$res)) {
                analyses$cases[[aname]]$devpars$res = analyses$devpars$res
            }
        }
    }
    cases[[name]]$analyses = analyses
}

do_one_case = function(name, case, ov_dir) {
    print(paste0("  Case: ", name))
    odir = file.path(ov_dir, name)
    dir.create(odir, showWarnings = FALSE)

    imm_ov = repOverlap(immdata$data, .method = case$method, .verbose = FALSE)
    vis_args = case$vis_args
    vis_args$.data = imm_ov
    p = do_call(vis, vis_args)

    ofig = file.path(odir, paste0(name, ".png"))
    png(ofig, width = case$devpars$width, height = case$devpars$height, res = case$devpars$res)
    print(p)
    dev.off()

    if (!is.null(case$analyses$cases) && length(case$analyses$cases) > 0) {
        for (aname in names(case$analyses$cases)) {
            print(paste0("    Analysis: ", aname))
            imm_ova = repOverlapAnalysis(imm_ov, .method = case$analyses$cases[[aname]]$method)
            avis_args = case$analyses$cases[[aname]]$vis_args
            avis_args$.data = imm_ova
            ap = do_call(vis, avis_args)
            if (aname == "DEFAULT") {
                aofig = file.path(odir, paste0(name, "-analysis.png"))
            } else {
                aofig = file.path(odir, paste0(name, "-", aname, "-analysis.png"))
            }
            png(aofig, width = case$analyses$cases[[aname]]$devpars$width, height = case$analyses$cases[[aname]]$devpars$height, res = case$analyses$cases[[aname]]$devpars$res)
            print(ap)
            dev.off()
        }
    }
}

ov_dir = file.path(outdir, "overlap")
dir.create(ov_dir, showWarnings = FALSE)

print("- Overlap")
for (name in names(cases)) {
    do_one_case(name, cases[[name]], ov_dir)
}
