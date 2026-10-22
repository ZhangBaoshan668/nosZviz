# ---------------------------------------------------------------------------
# alpha 
# ---------------------------------------------------------------------------

#  alpha 
alpha_table <- function(x, branches = NULL) {
  branches <- branches %||% intersect(x$branches, names(x$alpha))
  branches <- branches[vapply(branches, function(b) !is.null(x$alpha[[b]]), logical(1))]
  if (length(branches) == 0) nosz_abort(" alpha ")
  df <- do.call(rbind, lapply(branches, function(b) {
    d <- x$alpha[[b]]
    d$Branch <- b
    d
  }))
  df <- merge(df, x$design, by = "SampleID", all.x = TRUE)
  df$Branch <- factor(df$Branch, levels = branches)
  df$SampleID <- factor(df$SampleID, levels = x$design$SampleID)
  df
}

ALPHA_LABEL <- c(richness = "Richness", chao1 = "Chao1", ACE = "ACE",
                 shannon = "Shannon", simpson = "Simpson",
                 invsimpson = "Inverse Simpson")

#  -> 
alpha_lab <- function(i) {
  out <- ifelse(i %in% names(ALPHA_LABEL), ALPHA_LABEL[i], i)
  unname(out)
}

#  p 
alpha_test <- function(value, group, method = c("wilcox", "kruskal", "t.test")) {
  method <- match.arg(method)
  ok <- stats::complete.cases(value, group)
  value <- value[ok]; group <- factor(group[ok])
  if (nlevels(group) < 2 || length(value) < 3) return(list(p = NA_real_, name = "n.s."))
  p <- try({
    if (nlevels(group) == 2) {
      if (method == "t.test") stats::t.test(value ~ group)$p.value
      else stats::wilcox.test(value ~ group)$p.value
    } else {
      stats::kruskal.test(value ~ group)$p.value
    }
  }, silent = TRUE)
  if (inherits(p, "try-error")) p <- NA_real_
  nm <- if (nlevels(group) == 2) {
    if (method == "t.test") "Welch t-test" else "Wilcoxon test"
  } else "Kruskal-Wallis test"
  list(p = p, name = nm)
}

#' Alpha diversity boxplots
#'
#' @param x A `noszref` object.
#' @param index One or several alpha indices, e.g. `"shannon"` or
#'   `c("richness", "shannon")`.
#' @param branch Branches to show, default all with alpha results.
#' @param test `"wilcox"`, `"t.test"` or `"kruskal"`.
#' @param palette Palette name.
#' @param point Logical, show individual samples.
#' @param facet Logical, one panel per branch.
#' @return A ggplot object.
#' @export
plot_alpha_box <- function(x, index = "shannon", branch = NULL,
                           test = c("wilcox", "t.test", "kruskal"),
                           palette = "nosz", point = TRUE, facet = TRUE) {
  stopifnot(inherits(x, "noszref"))
  test <- match.arg(test)
  df <- alpha_table(x, branch)
  index <- intersect(index, colnames(df))
  if (length(index) == 0) nosz_abort(" alpha ",
                                     paste(setdiff(colnames(df), c("SampleID", "Group", "Branch")), collapse = ", "))
  long <- tidyr::pivot_longer(df, dplyr::all_of(index), names_to = "Index",
                              values_to = "Value")
  long$Index <- factor(long$Index, levels = index,
                       labels = ifelse(index %in% names(ALPHA_LABEL),
                                       ALPHA_LABEL[index], index))
  stats_df <- long %>%
    dplyr::group_by(Branch, Index) %>%
    dplyr::summarise(t = list(alpha_test(Value, Group, test)), .groups = "drop") %>%
    dplyr::mutate(p = vapply(t, function(z) z$p, numeric(1)),
                  lab = paste0(vapply(t, function(z) z$name, character(1)),
                               ": p = ", fmt_p(p)))
  p <- ggplot2::ggplot(long, ggplot2::aes(Group, Value, fill = Group)) +
    ggplot2::geom_boxplot(width = 0.58, outlier.shape = NA,
                          colour = "black", linewidth = 0.3, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = group_palette(df$Group, palette)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nosz(x_angle = 30, legend_position = "none")
  if (point) {
    p <- p + ggplot2::geom_jitter(width = 0.12, size = 1, shape = 21,
                                  fill = "white", colour = "black", stroke = 0.25)
  }
  p <- p + ggplot2::geom_text(data = stats_df,
                              ggplot2::aes(x = 1, y = Inf, label = lab),
                              inherit.aes = FALSE, vjust = 1.4, hjust = 0,
                              size = 1.9, colour = "grey20")
  if (length(index) > 1) {
    p <- p + ggplot2::facet_wrap(~Index, scales = "free_y", nrow = 1)
  }
  if (facet && length(unique(long$Branch)) > 1) {
    p <- p + ggplot2::facet_wrap(~Branch, scales = "free_y", nrow = 1)
  }
  p
}

#' Alpha diversity violin plots
#'
#' @inheritParams plot_alpha_box
#' @param add_box Logical, draw an inner boxplot.
#' @return A ggplot object.
#' @export
plot_alpha_violin <- function(x, index = "shannon", branch = NULL,
                              palette = "nosz", add_box = TRUE) {
  stopifnot(inherits(x, "noszref"))
  df <- alpha_table(x, branch)
  index <- intersect(index, colnames(df))
  if (length(index) == 0) nosz_abort(" alpha ")
  long <- tidyr::pivot_longer(df, dplyr::all_of(index), names_to = "Index",
                              values_to = "Value")
  long$Index <- factor(long$Index, levels = index,
                       labels = ifelse(index %in% names(ALPHA_LABEL),
                                       ALPHA_LABEL[index], index))
  p <- ggplot2::ggplot(long, ggplot2::aes(Group, Value, fill = Group)) +
    ggplot2::geom_violin(scale = "width", alpha = 0.6, colour = NA, width = 0.85) +
    ggplot2::geom_jitter(width = 0.1, size = 1, shape = 21, fill = "white",
                         colour = "black", stroke = 0.25) +
    ggplot2::scale_fill_manual(values = group_palette(df$Group, palette)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nosz(x_angle = 30, legend_position = "none")
  if (add_box) {
    p <- p + ggplot2::geom_boxplot(width = 0.14, fill = "white",
                                   outlier.shape = NA, linewidth = 0.25)
  }
  p + ggplot2::facet_wrap(~Index, scales = "free_y", nrow = 1)
}

#' Alpha diversity radar chart
#'
#' All indices are min-max normalised to 0--1 and group means are drawn as
#' closed polygons, which gives an immediate visual signature of each group's
#' diversity profile.
#'
#' @param x A `noszref` object.
#' @param branch Branch name; by default the branch with the largest OTU table.
#' @param indices Indices to include.
#' @param palette Palette name.
#' @param normalise `"minmax"` (per index) or `"none"`.
#' @return A ggplot object.
#' @export
plot_alpha_radar <- function(x, branch = NULL, indices = NULL,
                             palette = "nosz", normalise = c("minmax", "none")) {
  stopifnot(inherits(x, "noszref"))
  normalise <- match.arg(normalise)
  df <- alpha_table(x, branch)
  if (length(unique(df$Branch)) > 1) {
    df$Branch <- factor(df$Branch, levels = levels(df$Branch))
    df <- df[df$Branch == levels(df$Branch)[1], , drop = FALSE]
    nosz_info("radar  ", levels(droplevels(df$Branch))[1], "")
  }
  idx <- indices %||% intersect(c("richness", "chao1", "shannon", "simpson",
                                  "invsimpson"), colnames(df))
  idx <- intersect(idx, colnames(df))
  sub <- df[, c("SampleID", "Group", idx), drop = FALSE]
  if (normalise == "minmax") {
    for (i in idx) {
      r <- range(sub[[i]], na.rm = TRUE)
      sub[[i]] <- if (diff(r) == 0) 0.5 else (sub[[i]] - r[1]) / diff(r)
    }
  }
  agg <- sub %>% dplyr::group_by(Group) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(idx), ~ mean(.x, na.rm = TRUE)),
                     .groups = "drop")
  long <- tidyr::pivot_longer(agg, dplyr::all_of(idx), names_to = "Index",
                              values_to = "Value")
  long$Index <- factor(long$Index, levels = idx,
                       labels = ifelse(idx %in% names(ALPHA_LABEL), ALPHA_LABEL[idx], idx))
  # 
  base <- long %>% dplyr::filter(Index == levels(Index)[1]) %>% dplyr::mutate(Index = factor("_close", levels = c(levels(long$Index), "_close")))
  lv <- c(levels(long$Index), "_close")
  long$Index <- factor(as.character(long$Index), levels = lv)
  base$Index <- factor("_close", levels = lv)
  closed <- rbind(long, base)
  ggplot2::ggplot(closed, ggplot2::aes(Index, Value, group = Group,
                                       colour = Group, fill = Group)) +
    ggplot2::geom_polygon(alpha = 0.16, linewidth = 0.45) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::coord_polar() +
    ggplot2::scale_colour_manual(values = group_palette(closed$Group, palette)) +
    ggplot2::scale_fill_manual(values = group_palette(closed$Group, palette)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nosz(border = FALSE, legend_position = "right") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey88",
                                                            linewidth = 0.25))
}

#' Rarefaction curves
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param by_group Logical, draw group means with a confidence ribbon instead of
#'   one line per sample.
#' @param palette Palette name.
#' @param show_points Logical, overlay points on the curves.
#' @return A ggplot object.
#' @export
plot_rarefaction <- function(x, branch, by_group = TRUE, palette = "nosz",
                             show_points = FALSE) {
  stopifnot(inherits(x, "noszref"))
  d <- nosz_need(x, "rare", branch)
  long <- tidyr::pivot_longer(d, -Depth, names_to = "SampleID", values_to = "Richness")
  long <- merge(long, x$design, by = "SampleID", all.x = TRUE)
  long <- long[long$Depth > 0, , drop = FALSE]
  long$SampleID <- factor(long$SampleID, levels = x$design$SampleID)
  p <- ggplot2::ggplot(long, ggplot2::aes(Depth, Richness))
  if (by_group) {
    agg <- long %>% dplyr::group_by(Group, Depth) %>%
      dplyr::summarise(m = mean(Richness, na.rm = TRUE),
                       s = stats::sd(Richness, na.rm = TRUE),
                       .groups = "drop")
    agg$lo <- pmax(agg$m - agg$s, 0)
    agg$hi <- agg$m + agg$s
    p <- p +
      ggplot2::geom_ribbon(data = agg, ggplot2::aes(x = Depth, ymin = lo, ymax = hi,
                                                    fill = Group), alpha = 0.18,
                           colour = NA, inherit.aes = FALSE) +
      ggplot2::geom_line(data = agg, ggplot2::aes(x = Depth, y = m, colour = Group),
                         linewidth = 0.55, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = group_palette(long$Group, palette)) +
      ggplot2::scale_colour_manual(values = group_palette(long$Group, palette))
  } else {
    p <- p + ggplot2::geom_line(ggplot2::aes(colour = SampleID), linewidth = 0.4) +
      ggplot2::scale_colour_manual(values = sample_colours(x$design, palette))
  }
  if (show_points) p <- p + ggplot2::geom_point(size = 0.5, alpha = 0.5)
  p +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::labs(x = "Sequencing depth (reads)", y = "Observed richness",
                  title = paste0(branch, " rarefaction")) +
    theme_nosz(grid = TRUE, legend_position = "right")
}

#' Heatmap of alpha diversity indices
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param indices Indices to include.
#' @param scale `"column"` (z-score across samples per index) or `"none"`.
#' @param annotation Logical, annotate samples by group.
#' @param palette Diverging palette.
#' @return A `pheatmap` object.
#' @export
plot_alpha_heatmap <- function(x, branch, indices = NULL, scale = "column",
                               annotation = TRUE, palette = "RdYlBu") {
  stopifnot(inherits(x, "noszref"))
  d <- nosz_need(x, "alpha", branch)
  idx <- indices %||% intersect(c("richness", "chao1", "ACE", "shannon",
                                  "simpson", "invsimpson"), colnames(d))
  m <- df_to_matrix(d[, c("SampleID", idx)], 1)
  m <- t(m)
  m <- m[, order_samples(colnames(m), x$design), drop = FALSE]
  ann <- NULL
  if (annotation) {
    ann <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
    rownames(ann) <- ann$SampleID
    ann$SampleID <- NULL
    ann[] <- lapply(ann, as.character)
  }
  pheatmap::pheatmap(
    m, scale = if (scale == "none") "none" else "column",
    cluster_rows = nrow(m) > 2, cluster_cols = FALSE,
    annotation_col = ann, border_colour = "white", border_width = 0.8,
    display_numbers = TRUE, number_format = "%.2f", fontsize = 7,
    fontsize_number = 6,
    color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, palette))(101),
    angle_col = 45, main = paste0(branch, " alpha diversity (", scale, " scaled)")
  )
}

#' Scatter plot between two alpha diversity indices
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param index_x,index_y Indices for the two axes.
#' @param colour_by `"Group"` or `"SampleID"`.
#' @param smooth Logical, add a linear fit with confidence band.
#' @param label Logical, label samples.
#' @param palette Palette name.
#' @return A ggplot object.
#' @export
plot_alpha_scatter <- function(x, branch, index_x = "shannon",
                               index_y = "richness", colour_by = "Group",
                               smooth = TRUE, label = FALSE, palette = "nosz") {
  stopifnot(inherits(x, "noszref"))
  d <- nosz_need(x, "alpha", branch)
  if (!all(c(index_x, index_y) %in% colnames(d))) {
    nosz_abort("", paste(setdiff(colnames(d), "SampleID"), collapse = ", "))
  }
  d <- merge(d, x$design, by = "SampleID", all.x = TRUE)
  r <- suppressWarnings(stats::cor.test(d[[index_x]], d[[index_y]],
                                        method = "spearman"))
  lab <- paste0("Spearman r = ", formatC(unname(r$estimate), digits = 2, format = "f"),
                ", p = ", fmt_p(r$p.value))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data[[index_x]], .data[[index_y]],
                                       colour = .data[[colour_by]])) +
    ggplot2::geom_point(size = 1.6, alpha = 0.9) +
    ggplot2::labs(x = alpha_lab(index_x),
                  y = alpha_lab(index_y),
                  subtitle = lab, colour = colour_by) +
    theme_nosz(grid = TRUE, legend_position = "right")
  if (smooth) {
    p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                                  se = TRUE, alpha = 0.15, linewidth = 0.4)
  }
  if (label) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = SampleID),
                                      size = 1.8, seed = 1, min.segment.length = 0.2)
  }
  p + ggplot2::scale_colour_manual(values = group_palette(d[[colour_by]], palette))
}

#' Combined alpha diversity panel
#'
#' A six-index patchwork panel (richness, Chao1, ACE, Shannon, Simpson and
#' inverse Simpson) for one branch, ready for a supplementary figure.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param palette Palette name.
#' @param ncol Columns of the panel.
#' @return A patchwork object.
#' @export
plot_alpha_panel <- function(x, branch, palette = "nosz", ncol = 3) {
  stopifnot(inherits(x, "noszref"))
  d <- nosz_need(x, "alpha", branch)
  idx <- intersect(c("richness", "chao1", "ACE", "shannon", "simpson", "invsimpson"),
                   colnames(d))
  plist <- lapply(idx, function(i) {
    plot_alpha_box(x, index = i, branch = branch, palette = palette,
                   facet = FALSE) +
      ggplot2::labs(title = alpha_lab(i), y = alpha_lab(i)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"))
  })
  patchwork::wrap_plots(plist, ncol = ncol) +
    patchwork::plot_annotation(tag_levels = "A",
                               title = paste0(branch, " alpha diversity"))
}
