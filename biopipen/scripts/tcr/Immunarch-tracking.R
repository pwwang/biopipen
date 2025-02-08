log_info("")
log_info("# Clonotype tracking")
log_info("-----------------------------------")

log_info("Filling up cases ...")
if (is.null(trackings$subjects)) {
    trackings$subjects = c()
}
if (length(trackings$cases) == 0) {
    trackings$cases$DEFAULT = list(
        targets = trackings$targets,
        subject_col = trackings$subject_col,
        subjects = trackings$subjects
    )
} else {
    for (name in names(trackings$cases)) {
        if (is.null(trackings$cases[[name]]$targets)) {
            trackings$cases[[name]]$targets = trackings$targets
        }
        if (is.null(trackings$cases[[name]]$subject_col)) {
            trackings$cases[[name]]$subject_col = trackings$subject_col
        }
        if (is.null(trackings$cases[[name]]$subjects)) {
            trackings$cases[[name]]$subjects = trackings$subjects
        }
        if (is.null(trackings$cases[[name]]$subset)) {
            trackings$cases[[name]]$subset = trackings$subset
        }
    }
}

tracking_dir = file.path(outdir, "tracking")
dir.create(tracking_dir, showWarnings = FALSE)

run_tracking_case = function(casename) {
    case = trackings$cases[[casename]]

    if (!is.null(case$subset)) {
        d = immdata_from_expanded(filter_expanded_immdata(exdata, case$subset))
    } else {
        d = immdata
    }

    if (is.null(case$targets)) {
        # print(paste0("  ", casename, ", skip, no targets"))
        log_info("- Case: {casename}, skip, no targets")
    } else {
        # print(paste0("  ", casename))
        log_info("- Case: {casename}")
        allsubjects = d$meta %>% pull(case$subject_col) %>% unlist() %>% unique() %>% na.omit()
        if (is.null(case$subjects) || length(case$subjects) == 0) {
            subjects = allsubjects
        } else {
            subjects = intersect(case$subjects, allsubjects)
        }
        if (length(allsubjects) == 1) {
            stop(paste0("Cannot track clonotypes for only one subject: ", subjects))
        }
        samples = d$meta[d$meta[[case$subject_col]] %in% subjects, ] %>% pull(Sample) %>% unlist()
        if (is.numeric(case$targets)) {
            targets = do_call(rbind, lapply(samples, function(s) d$data[[s]])) %>%
                group_by(CDR3.aa) %>%
                summarise(Clones = sum(Clones)) %>%
                arrange(desc(Clones)) %>%
                slice_max(Clones, n=case$target) %>%
                pull(CDR3.aa)
        } else {
            targets = case$targets
        }
        if (case$subject_col == "Sample") {
            imm_tracking = trackClonotypes(d$data, targets, .col = "aa")
        } else {
            # Construct a data with names as subjects
            # For each subject, get the samples and merge the data
            # Then track
            newdata = list()
            for (subject in subjects) {
                subject_samples = d$meta[d$meta[[case$subject_col]] == subject, ] %>%
                    pull(Sample) %>%
                    na.omit()
                newdata[[subject]] = do_call(rbind, lapply(subject_samples, function(s) d$data[[s]]))
            }
            imm_tracking = trackClonotypes(newdata, targets, .col = "aa")
        }

        tracking_png = file.path(tracking_dir, paste0(slugify(casename), ".png"))
        png(tracking_png, res=100, height=1000, width=600 + 150 * length(subjects))
        print(vis(imm_tracking))
        dev.off()

        tracking_pdf = file.path(tracking_dir, paste0(slugify(casename), ".pdf"))
        pdf(tracking_pdf, height=10, width=6 + 1.5 * length(subjects))
        print(vis(imm_tracking))
        dev.off()

        add_report(
            list(
                kind = "descr",
                content = paste0(
                    "Clonotype tracking is a popular approach to monitor changes in the frequency of ",
                    "clonotypes of interest in vaccination and cancer immunology. ",
                    "For example, a researcher can track a clonotype across different time points ",
                    "in pre- and post-vaccination repertoires, or analyse the growth of ",
                    "malignant clonotypes in a tumor sample."
                )
            ),
            h1 = "Tracking of clonotypes"
        )

        add_report(
            list(
                src = tracking_png,
                download = tracking_pdf,
                name = if (casename == "DEFAULT") NULL else casename
            ),
            h1 = "Tracking of clonotypes",
            ui = "table_of_images"
        )
    }
}

sapply(names(trackings$cases), run_tracking_case)
