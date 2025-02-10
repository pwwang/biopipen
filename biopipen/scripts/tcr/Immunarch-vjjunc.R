log_info("")
log_info("# VJ Junction Circos Plots")
log_info("-----------------------------------")

# Already required by immunarch
library(circlize)

log_info("Filling up cases ...")
cases <- vj_juncs$cases
if (is.null(cases) || length(cases) == 0) {
    cases$DEFAULT <- list(
        by = vj_juncs$by,
        by_clones = vj_juncs$by_clones,
        subset = vj_juncs$subset,
        devpars = vj_juncs$devpars
    )
} else {
    for (name in names(cases)) {
        if (is.null(cases[[name]]$by)) {
            cases[[name]]$by <- vj_juncs$by
        }
        if (is.null(cases[[name]]$by_clones)) {
            cases[[name]]$by_clones <- vj_juncs$by_clones
        }
        if (is.null(cases[[name]]$subset)) {
            cases[[name]]$subset <- vj_juncs$subset
        }
        if (is.null(cases[[name]]$devpars)) {
            cases[[name]]$devpars <- vj_juncs$devpars
        }
        if (is.null(cases[[name]]$devpars$width)) {
            cases[[name]]$devpars$width <- vj_juncs$devpars$width
        }
        if (is.null(cases[[name]]$devpars$height)) {
            cases[[name]]$devpars$height <- vj_juncs$devpars$height
        }
        if (is.null(cases[[name]]$devpars$res)) {
            cases[[name]]$devpars$res <- vj_juncs$devpars$res
        }
    }
}

vjjunc_dir = file.path(outdir, "vj_junc")
dir.create(vjjunc_dir, showWarnings = FALSE)

do_one_case_vjjunc <- function(name, case) {
    log_info("Processing case: {name} ...")
    odir = file.path(vjjunc_dir, slugify(name))
    dir.create(odir, showWarnings = FALSE)

    if (!is.null(case$subset)) {
        d = filter_expanded_immdata(exdata, case$subset)
    } else {
        d = exdata
    }

    if (is.null(case$by) || length(case$by) == 0) {
        case$by <- "Sample"
    }

    by = trimws(strsplit(case$by, ",")[[1]])

    lapply(group_split(d, !!!syms(by)), function(gsd) {
        by_name <- gsd[1, by, drop = FALSE] %>% unlist() %>% paste0(collapse = "-")
        log_info("- Processing {by_name} ...")

        if (isTRUE(case$by_clones)) {
            gsd <- gsd %>% distinct(CDR3.aa, .keep_all = TRUE)
        }
        gsd <- gsd %>%
            group_by(V.name, J.name) %>%
            summarise(Size = n(), .groups = "drop") %>%
            filter(!is.na(V.name) & !is.na(J.name) & V.name != "None" & J.name != "None") %>%
            # if it's multiple chains, then split the chains
            separate_rows(V.name, J.name, sep = ";") %>%
            filter(!is.na(V.name) & !is.na(J.name) & V.name != "None" & J.name != "None") %>%
            group_by(V.name, J.name) %>%
            summarise(Size = sum(Size), .groups = "drop") %>%
            arrange(desc(Size), V.name, J.name)

        figfile <- file.path(odir, paste0(slugify(by_name), ".png"))
        png(figfile, width = case$devpars$width, height = case$devpars$height, res = case$devpars$res)
        circos.clear()
        tryCatch({
            chordDiagram(
                gsd,
                annotationTrack = c("grid", "axis"),
                preAllocateTracks = list(track.height = 0.25)
            )
        }, error = function(e) {
            log_warn("Error encountered: {e$message}, setting gap.after ...")
            circos.par(gap.after = c(rep(1, nrow(gsd) - 1), 5, rep(1, nrow(gsd) - 1), 5))
            chordDiagram(
                gsd,
                annotationTrack = c("grid", "axis"),
                preAllocateTracks = list(track.height = 0.25)
            )

        })
        circos.track(track.index = 1, panel.fun = function(x, y) {
            circos.text(
                CELL_META$xcenter,
                CELL_META$ylim[1],
                CELL_META$sector.index,
                cex = .8,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(-0.2, 0.5)
            )
        }, bg.border = NA) # here set bg.border to NA is important
        dev.off()

        # figfile_pdf <- file.path(odir, paste0(slugify(by_name), ".pdf"))
        # png(figfile_pdf, width = case$devpars$width / case$devpars$res, height = case$devpars$height / case$devpars$res)
        # circos.clear()
        # tryCatch({
        #     chordDiagram(
        #         gsd,
        #         annotationTrack = c("grid", "axis"),
        #         preAllocateTracks = list(track.height = 0.25)
        #     )
        # }, error = function(e) {
        #     log_warn("Error encountered: {e$message}, setting gap.after ...")
        #     circos.par(gap.after = c(rep(1, nrow(gsd) - 1), 5, rep(1, nrow(gsd) - 1), 5))
        #     chordDiagram(
        #         gsd,
        #         annotationTrack = c("grid", "axis"),
        #         preAllocateTracks = list(track.height = 0.25)
        #     )

        # })
        # circos.track(track.index = 1, panel.fun = function(x, y) {
        #     circos.text(
        #         CELL_META$xcenter,
        #         CELL_META$ylim[1],
        #         CELL_META$sector.index,
        #         cex = .8,
        #         facing = "clockwise",
        #         niceFacing = TRUE,
        #         adj = c(-0.2, 0.5)
        #     )
        # }, bg.border = NA) # here set bg.border to NA is important
        # dev.off()

        add_report(
            # list(src = figfile, name = by_name, download = figfile_pdf),
            list(src = figfile, name = by_name),
            h1 = "V-J Junction Circos Plots",
            h2 = ifelse(name == "DEFAULT", "#" , name),
            ui = "table_of_images"
        )

        NULL
    })
}

add_report(
    list(
        kind = "descr",
        content = "V-J usage plot displaying the frequency of various V-J junctions."
    ),
    h1 = "V-J Junction Circos Plots"
)

sapply(names(cases), function(name) do_one_case_vjjunc(name, cases[[name]]))
