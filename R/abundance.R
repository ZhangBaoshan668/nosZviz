# ---------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------

# 
BRANCH_COLS <- c(nosZI = "#1F78B4", nosZII = "#FF7F00", nosZIII = "#33A02C")

#' Branch-level abundance (RPKM)
#'
#' Visualizes the total abundance of nosZ-I / nosZ-II / nosZ-III in every
#' sample. Three layouts are available: per-sample bars, per-sample boxplots
#' and group summaries with error bars.
#'
#' @param x A `noszref` object.
#' @param style `"bar"` (default), `"box"` or `"summary"`.
#' @param log10 Logical, log10-transform the RPKM axis.
#' @param branches Branches to show, default all.
#' @param palette Palette name, see [nosz_pal()].
#' @param facet Logical, one panel per branch.
#' @param error `"sd"` or `"se"` for `style = "summary"`.
#' @param point Logical, overlay individual samples.
#' @return A ggplot object.
#' @export
plot_branch_rpkm <- function(x, style = c("bar", "box", "summary"),
                             log10 = TRUE, branches = NULL,
                             palette = "nosz", facet = TRUE, error = "sd",
                             point = TRUE) {
  stopifnot(inherits(x, "noszref"))
  style <- match.arg(style)
  branches <- branches %||% x$branches
  df <- do.call(rbind, lapply(branches, function(b) {
    d <- x$rpkm[[b]]
    if (is.null(d)) return(NULL)
    d$Branch <- b
    d
  }))
  if (is.null(df) || nrow(df) == 0) nosz_abort("No available branch RPKM table.")
  df <- merge(df, x$design, by = "SampleID", all.x = TRUE)
  df$SampleID <- factor(df$SampleID, levels = x$design$SampleID)
  df$Branch <- factor(df$Branch, levels = intersect(names(BRANCH_COLS), branches))
  if (log10) {
    df$Value <- log10(df$RPKM + 1)
    ylab <- expression(log[10]("RPKM + 1"))
  } else {
    df$Value <- df$RPKM
    ylab <- "RPKM"
  }

  if (style == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(SampleID, Value, fill = Group)) +
      ggplot2::geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
      ggplot2::labs(x = NULL, y = ylab, fill = "Group")
    if (point) {
      p <- p + ggplot2::geom_point(size = 0.9, colour = "black",
                                   position = ggplot2::position_jitter(width = 0.12))
    }
  } else if (style == "box") {
    p <- ggplot2::ggplot(df, ggplot2::aes(Group, Value, fill = Group)) +
      ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA,
                            colour = "black", linewidth = 0.3) +
      ggplot2::labs(x = NULL, y = ylab, fill = "Group")
    if (point) {
      p <- p + ggplot2::geom_jitter(width = 0.12, size = 1, shape = 21,
                                    colour = "black", stroke = 0.2)
    }
  } else {
    agg <- df %>%
      dplyr::group_by(Branch, Group) %>%
      dplyr::summarise(
        mean = mean(Value, na.rm = TRUE),
        sd = stats::sd(Value, na.rm = TRUE),
        n = dplyr::n(), .groups = "drop"
      ) %>%
      dplyr::mutate(err = if (error == "se") sd / sqrt(n) else sd)
    p <- ggplot2::ggplot(agg, ggplot2::aes(Group, mean, fill = Group)) +
      ggplot2::geom_col(width = 0.6, colour = "black", linewidth = 0.25) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - err, ymax = mean + err),
                             width = 0.18, linewidth = 0.3) +
      ggplot2::geom_point(data = df, ggplot2::aes(Group, Value),
                          inherit.aes = FALSE, shape = 21, fill = "white",
                          size = 1.2, stroke = 0.3,
                          position = ggplot2::position_jitter(width = 0.1)) +
      ggplot2::labs(x = NULL, y = ylab, fill = "Group")
  }

  p <- p + ggplot2::scale_fill_manual(values = group_palette(df$Group, palette)) +
    theme_nosz(x_angle = if (style == "bar") 45 else 0,
               legend_position = "top")
  if (facet && length(branches) > 1) {
    p <- p + ggplot2::facet_wrap(~Branch, scales = "free_y", nrow = 1)
  } else if (length(branches) == 1) {
    p <- p + ggplot2::labs(title = branches)
  }
  p
}

#' Branch composition of the nosZ community
#'
#' Shows how the total nosZ signal partitions among nosZ-I, nosZ-II and
#' nosZ-III, either as a 100% stacked bar per sample or as a donut chart of
#' group means.
#'
#' @param x A `noszref` object.
#' @param style `"stack"` or `"donut"`.
#' @param by `"sample"` or `"group"` (only for `style = "stack"`).
#' @param value `"rpkm"` (default) or `"count"` (read counts).
#' @param palette Palette name or `"branch"` for the built-in branch colours.
#' @return A ggplot object.
#' @export
plot_branch_composition <- function(x, style = c("stack", "donut"),
                                    by = c("sample", "group"),
                                    value = c("rpkm", "count"),
                                    palette = "branch") {
  stopifnot(inherits(x, "noszref"))
  style <- match.arg(style); by <- match.arg(by); value <- match.arg(value)
  src <- if (value == "rpkm") x$rpkm else x$counts
  df <- do.call(rbind, lapply(x$branches, function(b) {
    d <- src[[b]]
    if (is.null(d)) return(NULL)
    nm <- if ("RPKM" %in% names(d)) "RPKM" else "Count"
    data.frame(SampleID = d$SampleID, Branch = b, Value = d[[nm]],
               stringsAsFactors = FALSE)
  }))
  if (is.null(df)) nosz_abort("No available branch abundance table.")
  df <- merge(df, x$design, by = "SampleID", all.x = TRUE)
  df$Branch <- factor(df$Branch, levels = intersect(names(BRANCH_COLS), x$branches))

  cols <- if (identical(palette, "branch")) BRANCH_COLS else
    stats::setNames(nosz_pal(length(x$branches), palette), levels(df$Branch))

  if (style == "stack" && by == "sample") {
    p <- ggplot2::ggplot(df, ggplot2::aes(SampleID, Value, fill = Branch)) +
      ggplot2::geom_col(position = "fill", width = 0.72,
                        colour = "black", linewidth = 0.25) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::labs(x = NULL, y = "Relative abundance", fill = NULL) +
      ggplot2::scale_fill_manual(values = cols) +
      theme_nosz(legend_position = "top")
  } else if (style == "stack") {
    agg <- df %>%
      dplyr::group_by(Group, Branch) %>%
      dplyr::summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
    p <- ggplot2::ggplot(agg, ggplot2::aes(Group, Value, fill = Branch)) +
      ggplot2::geom_col(position = "fill", width = 0.6,
                        colour = "black", linewidth = 0.25) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::labs(x = NULL, y = "Relative abundance", fill = NULL) +
      ggplot2::scale_fill_manual(values = cols) +
      theme_nosz(legend_position = "top")
  } else {
    agg <- df %>%
      dplyr::group_by(Group, Branch) %>%
      dplyr::summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(Group) %>%
      dplyr::mutate(frac = Value / sum(Value),
                    ymax = cumsum(frac),
                    ymin = c(0, head(ymax, -1))) %>%
      dplyr::ungroup()
    p <- ggplot2::ggplot(agg) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 2, xmax = 3, ymin = ymin,
                                      ymax = ymax, fill = Branch),
                         colour = "white", linewidth = 0.6) +
      ggplot2::geom_text(ggplot2::aes(x = 3.25, y = (ymin + ymax) / 2,
                                      label = scales::percent(frac, accuracy = 0.1)),
                         hjust = 0, size = 2.1) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::xlim(0, 3.9) +
      ggplot2::facet_wrap(~Group) +
      ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
      ggplot2::scale_fill_manual(values = cols) +
      theme_nosz(border = FALSE, legend_position = "bottom") +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank())
  }
  p
}

#' Read partition Sankey / alluvial diagram
#'
#' Traces how raw sequencing reads are funnelled into the three nosZ branches
#' in each sample.
#'
#' @param x A `noszref` object.
#' @param palette Palette name or `"branch"`.
#' @param alpha Ribbon transparency.
#' @param label Logical, draw stratum labels.
#' @return A ggplot object.
#' @export
plot_reads_partition <- function(x, palette = "branch", alpha = 0.75,
                                 label = TRUE) {
  stopifnot(inherits(x, "noszref"))
  if (is.null(x$total_reads) || is.null(x$counts)) {
    nosz_abort(" total reads  reads ")
  }
  tr <- x$total_reads
  cnt <- do.call(rbind, lapply(x$branches, function(b) {
    d <- x$counts[[b]]
    if (is.null(d)) return(NULL)
    data.frame(SampleID = d$SampleID, Branch = b, Reads = d$Count,
               stringsAsFactors = FALSE)
  }))
  if (is.null(cnt)) nosz_abort("No branch reads table.")
  df <- merge(cnt, tr, by = "SampleID", all.x = TRUE)
  df <- merge(df, x$design, by = "SampleID", all.x = TRUE)
  unassigned <- df %>%
    dplyr::group_by(SampleID, Group, TotalReads) %>%
    dplyr::summarise(assigned = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Reads = pmax(TotalReads - assigned, 0), Branch = "unassigned")
  df$Branch <- as.character(df$Branch)
  df <- rbind(df[, c("SampleID", "Group", "Branch", "Reads")],
              unassigned[, c("SampleID", "Group", "Branch", "Reads")])
  df <- df[df$Reads > 0, , drop = FALSE]
  lv <- c("unassigned", intersect(names(BRANCH_COLS), x$branches))
  df$Branch <- factor(df$Branch, levels = lv)
  df$SampleID <- factor(df$SampleID, levels = x$design$SampleID)
  df <- df[order(df$SampleID, df$Branch), , drop = FALSE]

  cols <- if (identical(palette, "branch")) {
    c(unassigned = "grey80", BRANCH_COLS)
  } else {
    stats::setNames(nosz_pal(length(lv), palette), lv)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(y = Reads, axis1 = SampleID,
                                        axis2 = Branch)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = Branch), width = 1 / 5,
                              alpha = alpha, knot.pos = 0.35) +
    ggalluvial::geom_stratum(width = 1 / 5, fill = "white",
                             colour = "grey30", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = cols, name = "Branch") +
    ggplot2::scale_x_discrete(limits = c("Sample", "Branch"), expand = c(0.05, 0.05)) +
    ggplot2::labs(x = NULL, y = "Reads") +
    theme_nosz(legend_position = "top", grid = TRUE)
  if (label) {
    p <- p + ggplot2::geom_text(stat = ggalluvial::StatStratum,
                                ggplot2::aes(label = ggplot2::after_stat(stratum)),
                                size = 1.9, colour = "black")
  }
  p
}

#' Rank-abundance (Whittaker) curves
#'
#' @param x A `noszref` object.
#' @param branches Branches to draw, default all available.
#' @param by_sample Logical, one curve per sample (faceted by branch) instead
#'   of one mean curve per branch.
#' @param palette Palette name or `"branch"`.
#' @return A ggplot object.
#' @export
plot_rank_abundance <- function(x, branches = NULL, by_sample = FALSE,
                                palette = "branch") {
  stopifnot(inherits(x, "noszref"))
  branches <- branches %||% x$branches
  res <- lapply(branches, function(b) {
    m <- x$otutab[[b]]
    if (is.null(m)) return(NULL)
    if (by_sample) {
      long <- matrix_to_long(rel_abund(m), x$design)
      long <- long[long$Abundance > 0, , drop = FALSE]
      if (nrow(long) == 0) return(NULL)
      long <- long %>%
        dplyr::group_by(SampleID) %>%
        dplyr::arrange(dplyr::desc(Abundance), .by_group = TRUE) %>%
        dplyr::mutate(Rank = dplyr::row_number()) %>%
        dplyr::ungroup()
      long$Branch <- b
      long
    } else {
      v <- sort(rowMeans(rel_abund(m)), decreasing = TRUE)
      v <- v[v > 0]
      if (length(v) == 0) return(NULL)
      data.frame(Rank = seq_along(v), Abundance = v, Branch = b,
                 stringsAsFactors = FALSE)
    }
  })
  df <- do.call(rbind, res)
  if (is.null(df)) nosz_abort("No available OTU table.")
  df$Branch <- factor(df$Branch, levels = intersect(names(BRANCH_COLS), branches))
  cols <- if (identical(palette, "branch")) BRANCH_COLS else
    stats::setNames(nosz_pal(length(branches), palette), levels(df$Branch))

  p <- if (by_sample) {
    ggplot2::ggplot(df, ggplot2::aes(Rank, Abundance, group = SampleID,
                                     colour = Branch)) +
      ggplot2::geom_line(linewidth = 0.35, alpha = 0.8) +
      ggplot2::facet_wrap(~Branch, scales = "free")
  } else {
    ggplot2::ggplot(df, ggplot2::aes(Rank, Abundance, colour = Branch)) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::geom_point(size = 0.6, alpha = 0.6)
  }
  p + ggplot2::scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::annotation_logticks(sides = "bl", linewidth = 0.25) +
    ggplot2::labs(x = "OTU rank", y = "Relative abundance") +
    theme_nosz(grid = TRUE, legend_position = "top")
}

#' Heatmap of branch abundance across samples
#'
#' @param x A `noszref` object.
#' @param value `"rpkm"` or `"count"`.
#' @param scale `"row"`, `"column"` or `"none"`, passed to [pheatmap::pheatmap()].
#' @param show_numbers Logical, print values in cells.
#' @return A `pheatmap` object (print it to draw).
#' @export
plot_branch_heatmap <- function(x, value = c("rpkm", "count"), scale = "row",
                                show_numbers = TRUE) {
  stopifnot(inherits(x, "noszref"))
  value <- match.arg(value)
  src <- if (value == "rpkm") x$rpkm else x$counts
  col <- if (value == "rpkm") "RPKM" else "Count"
  df <- do.call(rbind, lapply(x$branches, function(b) {
    d <- src[[b]]
    if (is.null(d)) return(NULL)
    data.frame(SampleID = d$SampleID, Branch = b, Value = d[[col]],
               stringsAsFactors = FALSE)
  }))
  if (is.null(df)) nosz_abort("No available branch abundance table.")
  m <- long_to_matrix(df, "Branch", "SampleID", "Value")
  m <- m[, order_samples(colnames(m), x$design), drop = FALSE]
  ann <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
  rownames(ann) <- ann$SampleID
  ann$SampleID <- NULL
  ann[] <- lapply(ann, as.character)
  if (value == "rpkm") m <- log10(m + 1)
  pheatmap::pheatmap(
    m, scale = if (scale == "none") "none" else scale,
    cluster_rows = nrow(m) > 1, cluster_cols = FALSE,
    annotation_col = ann, border_colour = "white", border_width = 0.8,
    display_numbers = show_numbers, number_format = "%.1f",
    fontsize = 7, fontsize_number = 6,
    color = grDevices::colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(101),
    angle_col = 45,
    main = ifelse(value == "rpkm", "log10(RPKM+1), row-scaled", "Read counts")
  )
}

#' Violin + boxplot of branch abundance distribution
#'
#' @param x A `noszref` object.
#' @param value `"rpkm"` or `"count"`.
#' @param log10 Logical, log10 transform.
#' @param palette Palette name.
#' @return A ggplot object.
#' @export
plot_branch_violin <- function(x, value = c("rpkm", "count"), log10 = TRUE,
                               palette = "nosz") {
  stopifnot(inherits(x, "noszref"))
  value <- match.arg(value)
  src <- if (value == "rpkm") x$rpkm else x$counts
  col <- if (value == "rpkm") "RPKM" else "Count"
  df <- do.call(rbind, lapply(x$branches, function(b) {
    d <- src[[b]]
    if (is.null(d)) return(NULL)
    data.frame(SampleID = d$SampleID, Branch = b, Value = d[[col]],
               stringsAsFactors = FALSE)
  }))
  df <- merge(df, x$design, by = "SampleID", all.x = TRUE)
  df$Branch <- factor(df$Branch, levels = intersect(names(BRANCH_COLS), x$branches))
  ylab <- col
  if (log10) {
    df$Value <- log10(df$Value + 1)
    ylab <- paste0("log10(", col, " + 1)")
  }
  n_grp <- length(unique(df$Group))
  ggplot2::ggplot(df, ggplot2::aes(Branch, Value, fill = Branch)) +
    ggplot2::geom_violin(scale = "width", alpha = 0.55, colour = NA,
                         width = 0.8) +
    ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white",
                          linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(shape = Group), size = 1.1,
                        position = ggplot2::position_jitter(width = 0.07),
                        colour = "black") +
    ggplot2::scale_fill_manual(values = BRANCH_COLS) +
    ggplot2::labs(x = NULL, y = ylab) +
    theme_nosz(legend_position = if (n_grp > 1) "top" else "none")
}
