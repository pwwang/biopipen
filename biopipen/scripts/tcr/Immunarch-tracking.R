print("- Clonotype tracking")

trackings = {{ envs.trackings | r }}

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
    }
}

tracking_dir = file.path(outdir, "tracking")
dir.create(tracking_dir, showWarnings = FALSE)

run_tracking_case = function(casename) {
    case = trackings$cases[[casename]]

    if (is.null(case$targets)) {
        print(paste0("  ", casename, ", skip, no targets"))
    } else {
        print(paste0("  ", casename))
        allsubjects = immdata$meta %>% pull(case$subject_col) %>% unlist() %>% unique() %>% na.omit()
        if (is.null(case$subjects) || length(case$subjects) == 0) {
            subjects = allsubjects
        } else {
            subjects = intersect(case$subjects, allsubjects)
        }
        if (length(allsubjects) == 1) {
            stop(paste0("Cannot track clonotypes for only one subject: ", subjects))
        }
        samples = immdata$meta[immdata$meta[[case$subject_col]] %in% subjects, ] %>% pull(Sample) %>% unlist()
        if (is.numeric(case$targets)) {
            targets = do_call(rbind, lapply(samples, function(s) immdata$data[[s]])) %>%
                group_by(CDR3.aa) %>%
                summarise(Clones = sum(Clones)) %>%
                arrange(desc(Clones)) %>%
                slice_max(Clones, n=case$target) %>%
                pull(CDR3.aa)
        } else {
            targets = case$targets
        }
        if (case$subject_col == "Sample") {
            imm_tracking = trackClonotypes(immdata$data, targets, .col = "aa")
        } else {
            # Construct a data with names as subjects
            # For each subject, get the samples and merge the data
            # Then track
            newdata = list()
            for (subject in subjects) {
                subject_samples = immdata$meta[immdata$meta[[case$subject_col]] == subject, ] %>%
                    pull(Sample) %>%
                    na.omit()
                newdata[[subject]] = do_call(rbind, lapply(subject_samples, function(s) immdata$data[[s]]))
            }
            imm_tracking = trackClonotypes(newdata, targets, .col = "aa")
        }

        tracking_png = file.path(tracking_dir, paste0(casename, ".png"))
        png(tracking_png, res=100, height=1000, width=600 + 150 * length(subjects))
        print(vis(imm_tracking))
        dev.off()
    }
}

sapply(names(trackings$cases), run_tracking_case)
