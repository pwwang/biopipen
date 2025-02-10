# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

log_info("")
log_info("# Spectratyping analysis")
log_info("-----------------------------------")

# Fill up cases
log_info("Filling up cases ...")
if (is.null(spects$cases) || length(spects$cases) == 0) {
    spects$cases$DEFAULT = list(
        quant = spects$quant,
        col = spects$col,
        devpars = list(width = 1000, height = 1000, res = 150)
    )
} else {
    for (name in names(spects$cases)) {
        if (is.null(spects$cases[[name]]$quant)) {
            spects$cases[[name]]$quant = spects$quant
        }
        if (is.null(spects$cases[[name]]$col)) {
            spects$cases[[name]]$col = spects$col
        }
        if (is.null(spects$cases[[name]]$devpars)) {
            spects$cases[[name]]$devpars = list(width = 1000, height = 1000, res = 150)
        }
        if (is.null(spects$cases[[name]]$devpars$width)) {
            spects$cases[[name]]$devpars$width = 1000
        }
        if (is.null(spects$cases[[name]]$devpars$height)) {
            spects$cases[[name]]$devpars$height = 1000
        }
        if (is.null(spects$cases[[name]]$devpars$res)) {
            spects$cases[[name]]$devpars$res = 150
        }
        if (is.null(spects$cases[[name]]$subset)) {
            spects$cases[[name]]$subset = spects$subset
        }
    }
}

do_one_case_spectratyping = function(name, case, spect_dir) {
    # print(paste0("  Case: ", name))
    log_info("- Processing case: {name} ...")
    odir = file.path(spect_dir, slugify(name))
    dir.create(odir, showWarnings = FALSE)

    if (!is.null(case$subset)) {
        d = immdata_from_expanded(filter_expanded_immdata(exdata, case$subset))
    } else {
        d = immdata
    }

    for (sample in names(d$data)) {
        # print(paste0("    Sample: ", sample))
        log_info("  Sample: {sample} ...")
        spec_obj = spectratype(
            d$data[[sample]],
            .quant = case$quant,
            .col = case$col
        )
        spectfile = file.path(odir, paste0(slugify(sample), ".spect.png"))
        png(
            spectfile,
            res = case$devpars$res,
            width = case$devpars$width,
            height = case$devpars$height
        )
        print(vis(spec_obj))
        dev.off()

        spectfile_pdf = file.path(odir, paste0(slugify(sample), ".spect.pdf"))
        pdf(
            spectfile_pdf,
            width = case$devpars$width / case$devpars$res,
            height = case$devpars$height / case$devpars$res
        )
        print(vis(spec_obj))
        dev.off()

        add_report(
            list(src = spectfile, name = sample, download = spectfile_pdf),
            h1 = "Spectratyping",
            h2 = name,
            ui = "table_of_images"
        )
    }
}

add_report(
    list(
        kind = "descr",
        content = "Spectratype is a useful way to represent distributions of genes per sequence length."
    ),
    h1 = "Spectratyping"
)

spect_dir = file.path(outdir, "spectratyping")
dir.create(spect_dir, showWarnings = FALSE)

for (name in names(spects$cases)) {
    do_one_case_spectratyping(name, spects$cases[[name]], spect_dir)
}