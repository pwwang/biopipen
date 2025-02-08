# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

log_info("")
log_info("# K-mer analysis")
log_info("-----------------------------------")

# Fill up cases
log_info("Filling up cases ...")
cases = kmers$cases
if (is.null(cases) || length(cases) == 0) {
    cases$DEFAULT = list(
        k = kmers$k,
        head = kmers$head,
        vis_args = kmers$vis_args,
        devpars = kmers$devpars,
        profiles = kmers$profiles
    )
} else {
    for (name in names(cases)) {
        if (is.null(cases[[name]]$k)) {
            cases[[name]]$k = kmers$k
        }
        if (is.null(cases[[name]]$head)) {
            cases[[name]]$head = kmers$head
        }
        if (is.null(cases[[name]]$vis_args)) {
            cases[[name]]$vis_args = kmers$vis_args
        }
        if (is.null(cases[[name]]$devpars)) {
            cases[[name]]$devpars = kmers$devpars
        }
        if (is.null(cases[[name]]$devpars$width)) {
            cases[[name]]$devpars$width = kmers$devpars$width
        }
        if (is.null(cases[[name]]$devpars$height)) {
            cases[[name]]$devpars$height = kmers$devpars$height
        }
        if (is.null(cases[[name]]$devpars$res)) {
            cases[[name]]$devpars$res = kmers$devpars$res
        }
        if (is.null(cases[[name]]$profiles)) {
            cases[[name]]$profiles = kmers$profiles
        }
        if (is.null(cases[[name]]$subset)) {
            cases[[name]]$subset = kmers$subset
        }
    }
}

# profiles cases
for (name in names(cases)) {
    profiles = cases[[name]]$profiles
    if (is.null(profiles$cases) || length(profiles$cases) == 0) {
        profiles$cases$DEFAULT = list(
            method = profiles$method,
            vis_args = profiles$vis_args,
            devpars = profiles$devpars
        )
    } else {
        for (aname in names(profiles$cases)) {
            if (is.null(profiles$cases[[aname]]$method)) {
                profiles$cases[[aname]]$method = profiles$method
            }
            if (is.null(profiles$cases[[aname]]$vis_args)) {
                profiles$cases[[aname]]$vis_args = profiles$vis_args
            }
            if (is.null(profiles$cases[[aname]]$devpars)) {
                profiles$cases[[aname]]$devpars = profiles$devpars
            }
            if (is.null(profiles$cases[[aname]]$devpars$width)) {
                profiles$cases[[aname]]$devpars$width = profiles$devpars$width
            }
            if (is.null(profiles$cases[[aname]]$devpars$height)) {
                profiles$cases[[aname]]$devpars$height = profiles$devpars$height
            }
            if (is.null(profiles$cases[[aname]]$devpars$res)) {
                profiles$cases[[aname]]$devpars$res = profiles$devpars$res
            }
        }
    }
    cases[[name]]$profiles = profiles
}

do_one_case_kmer = function(name, case, kmer_dir) {
    # print(paste0("  Case: ", name))
    log_info("Processing case: {name} ...")
    odir = file.path(kmer_dir, slugify(name))
    dir.create(odir, showWarnings = FALSE)

    if (!is.null(case$subset)) {
        d = immdata_from_expanded(filter_expanded_immdata(exdata, case$subset))
    } else {
        d = immdata
    }

    imm_kmers = getKmers(d$data, case$k)
    vis_args = case$vis_args
    vis_args$.data = imm_kmers
    vis_args$.head = case$head
    p = do_call(vis, vis_args)

    ofig = file.path(odir, "Allsamples.png")
    png(ofig, width = case$devpars$width, height = case$devpars$height, res = case$devpars$res)
    print(p)
    dev.off()

    ofig_pdf = file.path(odir, "Allsamples.pdf")
    pdf(ofig_pdf, width = case$devpars$width / case$devpars$res, height = case$devpars$height / case$devpars$res)
    print(p)
    dev.off()

    add_report(
        list(
            kind = "descr",
            content = "K-mer sequence occurrences and motif analysis of CDR3 amino acid sequences"
        ),
        h1 = "Kmer and sequence motif analysis",
        h2 = ifelse(name == "DEFAULT", "#", name),
        h3 = "Kmer sequence occurrences"
    )

    add_report(
        list(kind = "image", src = ofig, download = ofig_pdf),
        h1 = "Kmer and sequence motif analysis",
        h2 = ifelse(name == "DEFAULT", "#", name),
        h3 = "Kmer sequence occurrences"
    )

    add_report(
        h1 = "Kmer and sequence motif analysis",
        h2 = ifelse(name == "DEFAULT", "#", name),
        h3 = "Motif analysis"
    )

    for (sample in names(d$data)) {
        # print(paste0("    Sample: ", sample))
        log_info("- Sample: {sample} ...")
        imm_kmer = getKmers(d$data[[sample]], case$k)

        if (!is.null(case$profiles$cases) && length(case$profiles$cases) > 0) {
            for (aname in names(case$profiles$cases)) {
                # print(paste0("    Profiling: ", aname))
                log_info("  Profiling: {aname} ...")
                imm_kmera = kmer_profile(imm_kmer, case$profiles$cases[[aname]]$method)
                avis_args = case$profiles$cases[[aname]]$vis_args
                avis_args$.data = imm_kmera
                ap = do_call(vis, avis_args)
                if (aname == "DEFAULT") {
                    aofig = file.path(odir, paste0(slugify(sample), "-profile.png"))
                } else {
                    aofig = file.path(odir, paste0(slugify(sample), "-", slugify(aname), "-profile.png"))
                }
                png(aofig, width = case$profiles$cases[[aname]]$devpars$width, height = case$profiles$cases[[aname]]$devpars$height, res = case$profiles$cases[[aname]]$devpars$res)
                print(ap)
                dev.off()

                aofig_pdf = gsub(".png$", ".pdf", aofig)
                pdf(aofig_pdf,
                    width = case$profiles$cases[[aname]]$devpars$width / case$profiles$cases[[aname]]$devpars$res,
                    height = case$profiles$cases[[aname]]$devpars$height / case$profiles$cases[[aname]]$devpars$res)
                print(ap)
                dev.off()

                add_report(
                    list(
                        src = aofig,
                        download = aofig_pdf,
                        name = paste0(sample, ifelse(aname == "DEFAULT", "", paste0(" - ", aname)))
                    ),
                    h1 = "Kmer and sequence motif analysis",
                    h2 = ifelse(name == "DEFAULT", "#", name),
                    h3 = "Motif analysis",
                    ui = "table_of_images"
                )
            }
        }
    }
}

add_report(
    list(
        kind = "descr",
        content = "Counting k-mer occurrences"
    ),
    h1 = "Kmer and sequence motif analysis"
)

kmer_dir = file.path(outdir, "kmer")
dir.create(kmer_dir, showWarnings = FALSE)

for (name in names(cases)) {
    do_one_case_kmer(name, cases[[name]], kmer_dir)
}
