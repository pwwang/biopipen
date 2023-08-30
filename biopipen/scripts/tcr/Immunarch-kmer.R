# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

kmers = {{ envs.kmers | r: todot="-" }}

# Fill up cases
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

do_one_case = function(name, case, kmer_dir) {
    print(paste0("  Case: ", name))
    odir = file.path(kmer_dir, name)
    dir.create(odir, showWarnings = FALSE)

    imm_kmers = getKmers(immdata$data, case$k)
    vis_args = case$vis_args
    vis_args$.data = imm_kmers
    vis_args$.head = case$head
    p = do_call(vis, vis_args)

    ofig = file.path(odir, "Allsamples.png")
    png(ofig, width = case$devpars$width, height = case$devpars$height, res = case$devpars$res)
    print(p)
    dev.off()

    for (sample in names(immdata$data)) {
        print(paste0("    Sample: ", sample))
        imm_kmer = getKmers(immdata$data[[sample]], case$k)

        if (!is.null(case$profiles$cases) && length(case$profiles$cases) > 0) {
            for (aname in names(case$profiles$cases)) {
                print(paste0("    Profiling: ", aname))
                imm_kmera = kmer_profile(imm_kmer, case$profiles$cases[[aname]]$method)
                avis_args = case$profiles$cases[[aname]]$vis_args
                avis_args$.data = imm_kmera
                ap = do_call(vis, avis_args)
                if (aname == "DEFAULT") {
                    aofig = file.path(odir, paste0(sample, "-profile.png"))
                } else {
                    aofig = file.path(odir, paste0(sample, "-", aname, "-profile.png"))
                }
                png(aofig, width = case$profiles$cases[[aname]]$devpars$width, height = case$profiles$cases[[aname]]$devpars$height, res = case$profiles$cases[[aname]]$devpars$res)
                print(ap)
                dev.off()
            }
        }
    }
}

kmer_dir = file.path(outdir, "kmer")
dir.create(kmer_dir, showWarnings = FALSE)

print("- K-mer analysis")
for (name in names(cases)) {
    do_one_case(name, cases[[name]], kmer_dir)
}
