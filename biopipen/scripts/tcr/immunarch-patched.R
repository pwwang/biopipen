library(immunarch)

vis.immunr_gini <- function(.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                            .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                            .legend = NA, .plot.type = "bar", ...) {
  # repDiversity(..., .method = "gini") generates a matrix
  .data = data.frame(Sample = rownames(.data), Value = .data[, 1])
  if (.plot.type == "bar") {
    vis_bar(
        .data = .data, .by = .by, .meta = .meta,
        .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
        .points = .points, .test = .test, .signif.label.size = .signif.label.size,
        .defgroupby = "Sample", .grouping.var = "Group",
        .labs = c(NA, "Gini coefficient"),
        .title = "Gini coefficient", .subtitle = "Sample diversity estimation using the Gini coefficient",
        .legend = .legend, .leg.title = NA
    )
  } else {
    vis_box(
        .data = .data, .by = .by, .meta = .meta, .test = .test,
        .points = .points, .signif.label.size = .signif.label.size,
        .defgroupby = "Sample", .grouping.var = "Group",
        .labs = c(NA, "Gini coefficient"),
        .title = "Gini coefficient", .subtitle = "Sample diversity estimation using the Gini coefficient",
        .legend = .legend, .leg.title = NA, .melt = FALSE
    )
  }
}

vis.immunr_div <- function(.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                            .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                            .legend = NA, .plot.type = "bar", ...) {
  # repDiversity(..., .method = "gini") generates a matrix
  if (.plot.type == "bar") {
    immunarch:::vis.immunr_div(.data = .data,.by = .by, .meta = .meta,
        .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
        .points = .points, .test = .test, .signif.label.size = .signif.label.size,
        .legend = .legend)
  } else {
    vis_box(
        .data = .data, .by = .by, .meta = .meta, .test = .test,
        .points = .points, .signif.label.size = .signif.label.size,
        .defgroupby = "Sample", .grouping.var = "Group",
        .labs = c(NA, "Effective number of clonoypes"),
        .title = "True diversity", .subtitle = "Sample diversity estimation using the true diversity index",
        .legend = NA, .leg.title = NA, .melt = FALSE
    )
  }
}

vis.immunr_chao1 <- function(.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                            .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                            .legend = NA, .plot.type = "bar", ...) {
  # repDiversity(..., .method = "gini") generates a matrix
  if (.plot.type == "bar") {
    immunarch:::vis.immunr_chao1(.data = .data,.by = .by, .meta = .meta,
        .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
        .points = .points, .test = .test, .signif.label.size = .signif.label.size,
        .legend = .legend)
  } else {
    .data <- data.frame(Sample = row.names(.data), Value = .data[, 1])
    vis_box(
        .data = .data, .by = .by, .meta = .meta, .test = .test,
        .points = .points, .signif.label.size = .signif.label.size,
        .defgroupby = "Sample", .grouping.var = "Group",
        .labs = c(NA, "Chao1"),
        .title = "Chao1", .subtitle = "Sample diversity estimation using Chao1",
        .legend = NA, .leg.title = NA, .melt = FALSE
    )
  }
}

vis.immunr_ginisimp <- function(.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                            .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                            .legend = NA, .plot.type = "bar", ...) {
  # repDiversity(..., .method = "gini") generates a matrix
  if (.plot.type == "bar") {
    immunarch:::vis.immunr_ginisimp(.data = .data,.by = .by, .meta = .meta,
        .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
        .points = .points, .test = .test, .signif.label.size = .signif.label.size,
        .legend = .legend)
  } else {
    vis_box(
        .data = .data, .by = .by, .meta = .meta, .test = .test,
        .points = .points, .signif.label.size = .signif.label.size,
        .defgroupby = "Sample", .grouping.var = "Group",
        .labs = c(NA, "Gini-Simpson index"),
        .title = "Gini-Simpson index", .subtitle = "Sample diversity estimation using the Gini-Simpson index",
        .legend = .legend, .leg.title = NA, .melt = FALSE
    )
  }
}

vis.immunr_invsimp <- function(.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                            .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                            .legend = NA, .plot.type = "bar", ...) {
  # repDiversity(..., .method = "gini") generates a matrix
  if (.plot.type == "bar") {
    immunarch:::vis.immunr_invsimp(.data = .data,.by = .by, .meta = .meta,
        .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
        .points = .points, .test = .test, .signif.label.size = .signif.label.size,
        .legend = .legend)
  } else {
    vis_box(
        .data = .data, .by = .by, .meta = .meta, .test = .test,
        .points = .points, .signif.label.size = .signif.label.size,
        .defgroupby = "Sample", .grouping.var = "Group",
        .labs = c(NA, "Inverse Simpson index"),
        .title = "Inverse Simpson index", .subtitle = "Sample diversity estimation using the inverse Simpson index",
        .legend = .legend, .leg.title = NA, .melt = FALSE
    )
  }
}

vis.immunr_dxx <- function(.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                            .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                            .legend = NA, .plot.type = "bar", ...) {
  # repDiversity(..., .method = "gini") generates a matrix
  if (.plot.type == "bar") {
    immunarch:::vis.immunr_dxx(.data = .data,.by = .by, .meta = .meta,
        .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
        .points = .points, .test = .test, .signif.label.size = .signif.label.size,
        .legend = .legend)
  } else {
    perc_value <- round(.data[1, 2][1])
    .data <- data.frame(Sample = row.names(.data), Value = .data[, 1])
    vis_box(
        .data = .data, .by = .by, .meta = .meta, .test = .test,
        .points = .points, .signif.label.size = .signif.label.size,
        .defgroupby = "Sample", .grouping.var = "Group",
        .labs = c(NA, paste0("D", perc_value)),
        .title = paste0("D", perc_value, " diversity index"), .subtitle = paste0("Number of clonotypes occupying the ", perc_value, "% of repertoires"),
        .legend = .legend, .leg.title = NA, .melt = FALSE
    )
  }
}
