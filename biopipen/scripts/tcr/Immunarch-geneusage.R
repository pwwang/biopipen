# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

log_info("")
log_info("# Gene usage analysis")
log_info("-----------------------------------")

# Fill up cases
log_info("Filling up cases ...")
cases = gene_usages$cases
if (is.null(cases) || length(cases) == 0) {
    cases$DEFAULT = list(
        top = gene_usages$top,
        norm = gene_usages$norm,
        by = gene_usages$by,
        vis_args = gene_usages$vis_args,
        devpars = gene_usages$devpars,
        analyses = gene_usages$analyses,
        subset = gene_usages$subset
    )
} else {
    for (name in names(cases)) {
        if (is.null(cases[[name]]$top)) {
            cases[[name]]$top = gene_usages$top
        }
        if (is.null(cases[[name]]$norm)) {
            cases[[name]]$norm = gene_usages$norm
        }
        if (is.null(cases[[name]]$by)) {
            cases[[name]]$by = gene_usages$by
        }
        if (is.null(cases[[name]]$vis_args)) {
            cases[[name]]$vis_args = gene_usages$vis_args
        }
        if (is.null(cases[[name]]$devpars)) {
            cases[[name]]$devpars = gene_usages$devpars
        }
        if (is.null(cases[[name]]$devpars$width)) {
            cases[[name]]$devpars$width = gene_usages$devpars$width
        }
        if (is.null(cases[[name]]$devpars$height)) {
            cases[[name]]$devpars$height = gene_usages$devpars$height
        }
        if (is.null(cases[[name]]$devpars$res)) {
            cases[[name]]$devpars$res = gene_usages$devpars$res
        }
        if (is.null(cases[[name]]$analyses)) {
            cases[[name]]$analyses = gene_usages$analyses
        }
        if (is.null(cases[[name]]$subset)) {
            cases[[name]]$subset = gene_usages$subset
        }
    }
}

# analyses cases
for (name in names(cases)) {
    if (is.character(cases[[name]]$by) && nchar(cases[[name]]$by) > 0) {
        cases[[name]]$by = trimws(unlist(strsplit(cases[[name]]$by, ",")))
    }
    analyses = cases[[name]]$analyses
    if (is.null(analyses$cases) || length(analyses$cases) == 0) {
        analyses$cases$DEFAULT = list(
            method = analyses$method,
            plot = analyses$plot,
            vis_args = analyses$vis_args,
            devpars = analyses$devpars
        )
    } else {
        for (aname in names(analyses$cases)) {
            if (is.null(analyses$cases[[aname]]$method)) {
                analyses$cases[[aname]]$method = analyses$method
            }
            if (is.null(analyses$cases[[aname]]$plot)) {
                analyses$cases[[aname]]$plot = analyses$plot
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

do_one_case_geneusage = function(name, case, gu_dir) {
    # print(paste0("  Case: ", name))
    log_info("Processing case: {name} ...")
    odir = file.path(gu_dir, name)
    dir.create(odir, showWarnings = FALSE)

    if (!is.null(case$subset)) {
        d = immdata_from_expanded(filter_expanded_immdata(exdata, case$subset))
    } else {
        d = immdata
    }

    imm_gu = geneUsage(d$data, .norm = case$norm)
    if (case$top > 0) {
        imm_gu = imm_gu %>%
            arrange(desc(rowSums(select(imm_gu, -"Names"), na.rm = TRUE))) %>%
            head(case$top)
    }

    vis_args = case$vis_args
    vis_args$.data = imm_gu
    if (!is.null(case$by)) {
        vis_args$.by = case$by
        vis_args$.meta = d$meta
    }
    p = do_call(vis, vis_args)

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
            kind = "table_image",
            src = ofig,
            download = ofig_pdf,
            descr = paste0(
                 "Distribution of known gene segments following the ",
                 '<a href="http://www.imgt.org/IMGTrepertoire/LocusGenes/" target="_blank">IMGT</a> ',
                 "nomenclature."
            )
        ),
        h1 = "Gene Usage",
        h2 = ifelse(name == "DEFAULT", "#", name)
    )

    if (!is.null(case$analyses$cases) && length(case$analyses$cases) > 0) {
        for (aname in names(case$analyses$cases)) {
            if (case$analyses$cases[[aname]]$method == "none") {
                next
            }
            # print(paste0("    Analysis: ", aname))
            log_info("- Processing analysis: {aname} ...")
            tryCatch({
                imm_gua = geneUsageAnalysis(imm_gu, .method = case$analyses$cases[[aname]]$method)
            }, error = function(e) {
                if (grepl("cmdscale", e)) {
                    stop(paste0(
                        "Too few samples (", n_samples, ") for gene usage analyses.\n",
                        "You can use a different method or set it to `none` to skip it.\n",
                        "(Orginal error: ", e, ")"
                    ))
                } else {
                    stop(e)
                }
            })
            avis_args = case$analyses$cases[[aname]]$vis_args
            avis_args$.data = imm_gua
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
                h1 = "Gene Usage",
                h2 = ifelse(name == "DEFAULT", "#", name),
                h3 = "Gene Usage Analysis",
                ui = "table_of_images"
            )
        }
    }
}

gu_dir = file.path(outdir, "gene_usage")
dir.create(gu_dir, showWarnings = FALSE)

for (name in names(cases)) {
    do_one_case_geneusage(name, cases[[name]], gu_dir)
}
