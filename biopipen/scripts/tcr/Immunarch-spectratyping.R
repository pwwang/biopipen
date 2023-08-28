# loaded variables
# immfile, outdir, mutaters, immdata, n_samples

spects = {{ envs.spects | r }}

# Fill up cases
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
    }
}

do_one_case = function(name, case, spect_dir) {
    print(paste0("  Case: ", name))
    odir = file.path(spect_dir, name)
    dir.create(odir, showWarnings = FALSE)

    for (sample in names(immdata$data)) {
        print(paste0("    Sample: ", sample))
        spec_obj = spectratype(
            immdata$data[[sample]],
            .quant = case$quant,
            .col = case$col
        )
        png(
            file.path(odir, paste0(sample, ".png")),
            res = case$devpars$res,
            width = case$devpars$width,
            height = case$devpars$height
        )
        print(vis(spec_obj))
        dev.off()
    }
}

print("- Spectratyping")
spect_dir = file.path(outdir, "spectratyping")
dir.create(spect_dir, showWarnings = FALSE)

for (name in names(spects$cases)) {
    do_one_case(name, spects$cases[[name]], spect_dir)
}