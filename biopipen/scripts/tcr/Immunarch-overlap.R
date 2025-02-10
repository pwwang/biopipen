# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

log_info("")
log_info("# Overlap analysis")
log_info("-----------------------------------")

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
        if (is.null(cases[[name]]$subset)) {
            cases[[name]]$subset = overlaps$subset
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

get_method_descr <- function(method) {
    descr <- switch(method,
        public = paste0(
            "number of public (shared) clonotypes, ",
            "a classic measure of overlap similarity"
        ),
        overlap = paste0(
            "overlap coefficient, a normalised measure of overlap similarity. ",
            "It is defined as the size of the intersection divided by the smaller of ",
            "the size of the two sets."
        ),
        jaccard = paste0(
            "Jaccard index, measures the similarity between finite sample sets, ",
            "and is defined as the size of the intersection divided by the size of ",
            "the union of the sample sets."
        ),
        tversky = paste0(
            "Tversky index, an asymmetric similarity measure on sets that compares ",
            "a variant to a prototype. ",
            "If using default arguments, it’s similar to Dice’s coefficient."
        ),
        cosine = "cosine similarity, a measure of similarity between two non-zero vectors",
        morisita = paste0(
            "Morisita's overlap index, a statistical measure of dispersion of ",
            "individuals in a population. ",
            "It is used to compare overlap among samples."
        )
    )

    if (!is.null(descr)) {
        return(descr)
    }

    return(paste0(
        "incremental overlap, ",
        "overlaps of the N most abundant clonotypes with incrementally growing N"
    ))
}

do_one_case_overlap = function(name, case, ov_dir) {
    # print(paste0("  Case: ", name))
    log_info("Processing case: {name} ...")
    odir = file.path(ov_dir, name)
    dir.create(odir, showWarnings = FALSE)

    if (!is.null(case$subset)) {
        d = immdata_from_expanded(filter_expanded_immdata(exdata, case$subset))
    } else {
        d = immdata
    }

    imm_ov = repOverlap(d$data, .method = case$method, .verbose = FALSE)
    vis_args = case$vis_args
    vis_args$.data = imm_ov
    p = do_call(vis, vis_args)

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
            kind = "table_image",
            src = ofig,
            download = ofig_pdf,
            descr = paste0(
                "Repertoire overlap is the most common approach to measure repertoire similarity, ",
                "using method <code>", case$method, "</code>, ",
                get_method_descr(case$method)
            )
        ),
        h1 = "Repertoire Overlaps",
        h2 = ifelse(name == "DEFAULT", "#", name)
    )

    if (!is.null(case$analyses$cases) && length(case$analyses$cases) > 0) {
        for (aname in names(case$analyses$cases)) {
            if (case$analyses$cases[[aname]]$method == "none") next

            # print(paste0("    Analysis: ", aname))
            log_info("- Processing analysis: {aname} ...")
            k = min(n_samples - 1, 2)
            tryCatch({
                imm_ova = repOverlapAnalysis(
                    imm_ov, .k = k, .perp = 1e-5, .method = case$analyses$cases[[aname]]$method
                )
            }, error = function(e) {
                if (grepl("cmdscale", e)) {
                    stop(paste0(
                        "Too few samples (", n_samples, ") for overlap analyses.\n",
                        "You can use a different method or set it to `none` to skip it.\n",
                        "(Orginal error: ", e, ")"
                    ))
                } else {
                    stop(e)
                }
            })
            avis_args = case$analyses$cases[[aname]]$vis_args
            avis_args$.data = imm_ova
            ap = do_call(vis, avis_args)
            if (aname == "DEFAULT") {
                aofig = file.path(odir, paste0(name, "-analysis.png"))
                aofig_pdf = file.path(odir, paste0(name, "-analysis.pdf"))
            } else {
                aofig = file.path(odir, paste0(name, "-", aname, "-analysis.png"))
                aofig_pdf = file.path(odir, paste0(name, "-", aname, "-analysis.pdf"))
            }
            png(aofig, width = case$analyses$cases[[aname]]$devpars$width, height = case$analyses$cases[[aname]]$devpars$height, res = case$analyses$cases[[aname]]$devpars$res)
            print(ap)
            dev.off()

            pdf(aofig_pdf,
                width = case$analyses$cases[[aname]]$devpars$width / case$analyses$cases[[aname]]$devpars$res,
                height = case$analyses$cases[[aname]]$devpars$height / case$analyses$cases[[aname]]$devpars$res)
            print(ap)
            dev.off()

            add_report(
                list(src = aofig, name = aname, download = aofig_pdf),
                h1 = "Repertoire Overlaps",
                h2 = ifelse(name == "DEFAULT", "#", name),
                h3 = "Repertoire Overlap Analysis",
                ui = "table_of_images"
            )

        }
    }
}

ov_dir = file.path(outdir, "overlap")
dir.create(ov_dir, showWarnings = FALSE)

# print("- Overlap")
for (name in names(cases)) {
    do_one_case_overlap(name, cases[[name]], ov_dir)
}
