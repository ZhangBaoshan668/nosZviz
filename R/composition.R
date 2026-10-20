# ---------------------------------------------------------------------------
# 模块二：分类组成
# ---------------------------------------------------------------------------

# 内部分支图：relative abundance 长表
tax_rel_long <- function(x, branch, level, top_n = 10) {
  m <- nosz_tax_matrix(x, branch, level)
  r <- rel_abund(m)
  r <- collapse_taxa(r, top_n = top_n)
  long <- matrix_to_long(r, x$design)
  long$Taxon <- factor(as.character(long$Taxon), levels = rownames(r))
  long$Branch <- branch
  long
}

# 分类单元配色：Other 固定灰色
tax_colours <- function(taxa, palette = "Paired") {
  taxa <- as.character(taxa)
  cols <- nosz_pal(length(taxa), palette)
  names(cols) <- taxa
  if ("Other" %in% taxa) cols["Other"] <- "grey80"
  cols
}

#' Taxonomic composition stacked barplot
#'
#' The classic community composition plot, with support for sample-level and
#' group-mean layouts, automatic merging of low-abundance taxa into `Other`,
#' and multiple colour palettes. When several branches are requested, one
#' panel per branch is returned as a patchwork object.
#'
#' @param x A `noszref` object.
#' @param level Taxonomic level.
#' @param branch One branch, or `NULL` for all available branches.
#' @param style `"sample"` (default) or `"group"`.
#' @param top_n Number of most abundant taxa kept, the rest become `Other`.
#' @param palette Palette name, see [nosz_pal()].
#' @param x_angle Rotation of x labels (group style).
#' @return A ggplot or patchwork object.
#' @export
plot_tax_stack <- function(x, level = "phylum", branch = NULL,
                           style = c("sample", "group"), top_n = 10,
                           palette = "Paired", x_angle = 45) {
  stopifnot(inherits(x, "noszref"))
  style <- match.arg(style)
  branches <- branch %||% intersect(x$branches, names(x$tax))
  branches <- branches[vapply(branches, function(b) !is.null(x$tax[[b]][[level]]), logical(1))]
  if (length(branches) == 0) nosz_abort("没有分支包含 ", level, " 层级。")

  one <- function(b) {
    long <- tax_rel_long(x, b, level, top_n)
    if (style == "group") {
      long <- long %>%
        dplyr::group_by(Group, Taxon) %>%
        dplyr::summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
      long$x <- long$Group
    } else {
      long$x <- long$SampleID
    }
    cols <- tax_colours(levels(droplevels(long$Taxon)), palette)
    p <- ggplot2::ggplot(long, ggplot2::aes(x, Abundance, fill = Taxon)) +
      ggplot2::geom_col(width = 0.72, colour = "black", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = cols, name = RANK_LABEL[[level]]) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                  expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = "Relative abundance",
                    title = paste0(b, " - ", RANK_LABEL[[level]])) +
      theme_nosz(x_angle = if (style == "group") x_angle else 45,
                 legend_position = "right")
    if (style == "sample") {
      p <- p + ggplot2::scale_x_discrete(limits = x$design$SampleID)
    }
    p
  }

  plist <- lapply(branches, one)
  if (length(plist) == 1) return(plist[[1]])
  patchwork::wrap_plots(plist, ncol = 1, guides = "collect") &
    ggplot2::theme(legend.position = "right")
}

#' Heatmap of taxon abundance
#'
#' Row-scaled heatmap of the most abundant taxa with sample annotation.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param top_n Number of taxa to display.
#' @param scale `"row"`, `"column"` or `"none"`.
#' @param value `"rpkm"` or `"count"`.
#' @param cluster Logical, cluster rows.
#' @param annotation Logical, add sample group annotation bar.
#' @param palette Colour palette of the heatmap body.
#' @return A `pheatmap` object.
#' @export
plot_tax_heatmap <- function(x, branch, level = "genus", top_n = 30,
                             scale = "row", value = c("rpkm", "count"),
                             cluster = TRUE, annotation = TRUE,
                             palette = "RdYlBu") {
  stopifnot(inherits(x, "noszref"))
  value <- match.arg(value)
  m <- nosz_tax_matrix(x, branch, level, value = value)
  if (nrow(m) > top_n) {
    m <- m[order(rowMeans(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  }
  m <- m[, order_samples(colnames(m), x$design), drop = FALSE]
  if (value == "rpkm") m <- log10(m + 1)
  ann <- NULL
  if (annotation) {
    ann <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
    rownames(ann) <- ann$SampleID
    ann$SampleID <- NULL
    ann[] <- lapply(ann, as.character)
  }
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, palette))(101)
  pheatmap::pheatmap(
    m, scale = if (scale == "none") "none" else scale,
    cluster_rows = cluster && nrow(m) > 2, cluster_cols = FALSE,
    annotation_col = ann, color = pal, border_colour = "grey90",
    border_width = 0.4, fontsize = 6.5, fontsize_row = 6,
    angle_col = 45, treeheight_row = 12,
    main = paste0(branch, " | ", RANK_LABEL[[level]],
                  if (value == "rpkm") " | log10(RPKM+1)" else " | counts")
  )
}

#' Bubble plot of taxon abundance
#'
#' Encodes abundance as bubble size and row-scaled abundance as colour,
#' which is far more readable than a heatmap when there are many zero cells.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param top_n Number of taxa.
#' @param colour_by `"zscore"` or `"abundance"`.
#' @param palette Continuous palette name, see [nosz_pal()] or `"viridis"`.
#' @return A ggplot object.
#' @export
plot_tax_bubble <- function(x, branch, level = "genus", top_n = 20,
                            colour_by = c("zscore", "abundance"),
                            palette = "RdYlBu") {
  stopifnot(inherits(x, "noszref"))
  colour_by <- match.arg(colour_by)
  m <- nosz_tax_matrix(x, branch, level)
  if (nrow(m) > top_n) {
    m <- m[order(rowMeans(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  }
  m <- rel_abund(m) * 100
  z <- zscore_rows(m)
  long <- matrix_to_long(m, x$design)
  long$Z <- as.numeric(z)
  long$Taxon <- factor(long$Taxon, levels = rev(rownames(m)))
  long$SampleID <- factor(long$SampleID, levels = intersect(x$design$SampleID, colnames(m)))
  long$Value <- if (colour_by == "zscore") long$Z else long$Abundance
  pal <- if (identical(palette, "viridis")) {
    viridisLite::viridis(101)
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, palette))(101)
  }
  ggplot2::ggplot(long, ggplot2::aes(SampleID, Taxon)) +
    ggplot2::geom_point(ggplot2::aes(size = Abundance, fill = Value),
                        shape = 21, colour = "grey30", stroke = 0.2) +
    ggplot2::scale_size_area(max_size = 4.2, name = "Relative\nabundance (%)") +
    ggplot2::scale_fill_gradientn(
      colours = pal, name = if (colour_by == "zscore") "Row z-score" else "Abundance (%)"
    ) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0(branch, " - ", RANK_LABEL[[level]])) +
    theme_nosz(x_angle = 45, legend_position = "right", grid = TRUE)
}

#' Alluvial / streamgraph of taxonomic composition
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param top_n Number of taxa kept.
#' @param palette Palette name.
#' @param alpha Ribbon transparency.
#' @param style `"alluvium"` (flow between samples) or `"area"` (stacked area).
#' @return A ggplot object.
#' @export
plot_tax_area <- function(x, branch, level = "phylum", top_n = 8,
                          palette = "Paired", alpha = 0.85,
                          style = c("alluvium", "area")) {
  stopifnot(inherits(x, "noszref"))
  style <- match.arg(style)
  long <- tax_rel_long(x, branch, level, top_n)
  long$SampleID <- factor(long$SampleID, levels = intersect(x$design$SampleID,
                                                            unique(long$SampleID)))
  long <- long[order(long$SampleID, long$Taxon), , drop = FALSE]
  cols <- tax_colours(levels(droplevels(long$Taxon)), palette)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = SampleID, y = Abundance,
                                          alluvium = Taxon, stratum = Taxon,
                                          fill = Taxon)) +
    ggplot2::scale_fill_manual(values = cols, name = RANK_LABEL[[level]]) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Relative abundance",
                  title = paste0(branch, " - ", RANK_LABEL[[level]])) +
    theme_nosz(x_angle = 45, legend_position = "right")
  if (style == "alluvium") {
    p <- p + ggalluvial::geom_flow(alpha = alpha, colour = "white",
                                   linewidth = 0.15, curve_type = "sine") +
      ggalluvial::geom_stratum(width = 0.55, alpha = alpha,
                               colour = "grey25", linewidth = 0.15)
  } else {
    p <- p + ggplot2::geom_area(ggplot2::aes(group = Taxon),
                                colour = "white", linewidth = 0.15, alpha = alpha)
  }
  p
}

#' Nested donut (sunburst) chart of taxonomic composition
#'
#' Two concentric rings: the inner ring shows phylum-level fractions and the
#' outer ring the genus-level breakdown within each phylum, averaged per group.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param inner Inner taxonomic level, default `"phylum"`.
#' @param outer Outer taxonomic level, default `"genus"`.
#' @param top_outer Number of outer-ring taxa kept per inner group.
#' @param palette Palette name.
#' @param by `"group"` or `"sample"` (one ring per sample).
#' @return A ggplot object.
#' @export
plot_tax_sunburst <- function(x, branch, inner = "phylum", outer = "genus",
                              top_outer = 4, palette = "Paired",
                              by = c("group", "sample")) {
  stopifnot(inherits(x, "noszref"))
  by <- match.arg(by)
  mi <- nosz_tax_matrix(x, branch, inner)
  mo <- nosz_tax_matrix(x, branch, outer)
  # 将 outer 归属到 inner：用 BLAST 长表或按命名无法直接映射，
  # 这里按样本内 outer 排序做近似展示，并用 inner 总量控制比例
  smp <- intersect(colnames(mi), colnames(mo))
  if (length(smp) == 0) nosz_abort("内外层级样本不匹配。")
  unit <- if (by == "group") x$design$Group[match(smp, x$design$SampleID)] else smp
  mi <- rel_abund(mi[, smp, drop = FALSE])
  mo <- rel_abund(mo[, smp, drop = FALSE])
  grp <- factor(unit, levels = unique(unit))

  # inner 汇总
  di <- data.frame(
    Unit = rep(as.character(grp), each = nrow(mi)),
    Taxon = rep(rownames(mi), times = length(grp)),
    V = as.numeric(mi)
  )
  di <- di %>% dplyr::group_by(Unit, Taxon) %>% dplyr::summarise(V = mean(V), .groups = "drop")
  di <- di %>% dplyr::group_by(Unit) %>% dplyr::mutate(f = V / sum(V)) %>%
    dplyr::arrange(dplyr::desc(f), .by_group = TRUE) %>% dplyr::ungroup()
  di <- di %>% dplyr::group_by(Unit) %>% dplyr::mutate(ymax = cumsum(f), ymin = dplyr::lag(ymax, default = 0)) %>% dplyr::ungroup()

  do <- data.frame(
    Unit = rep(as.character(grp), each = nrow(mo)),
    Taxon = rep(rownames(mo), times = length(grp)),
    V = as.numeric(mo)
  )
  do <- do %>% dplyr::group_by(Unit, Taxon) %>% dplyr::summarise(V = mean(V), .groups = "drop")
  do <- do %>% dplyr::group_by(Unit) %>%
    dplyr::mutate(rank = dplyr::min_rank(dplyr::desc(V)),
                  Taxon = ifelse(rank > top_outer, "Other", Taxon)) %>%
    dplyr::group_by(Unit, Taxon) %>% dplyr::summarise(V = sum(V), .groups = "drop") %>%
    dplyr::group_by(Unit) %>% dplyr::mutate(f = V / sum(V)) %>%
    dplyr::arrange(dplyr::desc(f), .by_group = TRUE) %>% dplyr::ungroup()
  do <- do %>% dplyr::group_by(Unit) %>% dplyr::mutate(ymax = cumsum(f), ymin = dplyr::lag(ymax, default = 0)) %>% dplyr::ungroup()

  # 外环颜色继承内环：按 inner 名称匹配同名列，匹配不到用统一调色板
  inner_cols <- tax_colours(unique(di$Taxon), palette)
  outer_taxa <- unique(do$Taxon)
  outer_cols <- tax_colours(outer_taxa, "Set3")
  shared <- intersect(outer_taxa, names(inner_cols))
  outer_cols[shared] <- inner_cols[shared]

  di$Taxon <- factor(di$Taxon, levels = unique(di$Taxon))
  do$Taxon <- factor(do$Taxon, levels = unique(do$Taxon))

  ggplot2::ggplot() +
    ggplot2::geom_rect(data = di, ggplot2::aes(xmin = 1.6, xmax = 2.4,
                                               ymin = ymin, ymax = ymax,
                                               fill = Taxon),
                       colour = "white", linewidth = 0.4) +
    ggplot2::geom_rect(data = do, ggplot2::aes(xmin = 2.45, xmax = 3.3,
                                               ymin = ymin, ymax = ymax,
                                               fill = Taxon),
                       colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(data = di, ggplot2::aes(x = 2.0, y = (ymin + ymax) / 2,
                                               label = ifelse(f > 0.06, as.character(Taxon), "")),
                       size = 1.7, colour = "white", fontface = "bold") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(~Unit) +
    ggplot2::xlim(0, 3.9) +
    ggplot2::scale_fill_manual(values = c(inner_cols, outer_cols),
                               breaks = levels(di$Taxon), name = RANK_LABEL[[inner]]) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0(branch, ": ", RANK_LABEL[[inner]], " (inner) and ",
                                 RANK_LABEL[[outer]], " (outer)")) +
    theme_nosz(border = FALSE, legend_position = "right") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
}

#' Ternary plot of nosZ branch contribution per taxon
#'
#' Each vertex is one nosZ branch; every point is a taxon whose position
#' reflects its relative abundance among the three branches.
#'
#' @param x A `noszref` object.
#' @param level Taxonomic level.
#' @param top_n Number of taxa (by total abundance across branches).
#' @param colour_by Colour points by the dominant phylum (`"phylum"`) or by
#'   total abundance (`"abundance"`).
#' @param palette Palette name.
#' @param label Logical, label taxa with `ggrepel`.
#' @return A ggplot object.
#' @export
plot_tax_ternary <- function(x, level = "genus", top_n = 25,
                             colour_by = c("phylum", "abundance"),
                             palette = "Paired", label = TRUE) {
  stopifnot(inherits(x, "noszref"))
  colour_by <- match.arg(colour_by)
  need <- c("nosZI", "nosZII", "nosZIII")
  if (!all(need %in% names(x$tax))) {
    nosz_abort("三元图需要 nosZI、nosZII、nosZIII 三个分支的分类表。当前：",
               paste(names(x$tax), collapse = ", "))
  }
  mats <- lapply(need, function(b) {
    m <- try(nosz_tax_matrix(x, b, level), silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    rowSums(rel_abund(m))
  })
  if (any(vapply(mats, is.null, logical(1)))) nosz_abort("某些分支缺少该层级分类表。")
  taxa <- Reduce(union, lapply(mats, names))
  tab <- data.frame(Taxon = taxa, stringsAsFactors = FALSE)
  for (i in seq_along(need)) {
    v <- mats[[i]][taxa]
    v[is.na(v)] <- 0
    tab[[need[i]]] <- as.numeric(v)
  }
  tot <- rowSums(tab[, need])
  tab <- tab[tot > 0, , drop = FALSE]
  if (nrow(tab) > top_n) {
    tab <- tab[order(rowSums(tab[, need]), decreasing = TRUE)[seq_len(top_n)], ]
  }
  tab$frac <- tab[, need] / rowSums(tab[, need])
  # 三角坐标：A 顶点 (0,0)、B 顶点 (1,0)、C 顶点 (0.5, sqrt(3)/2)
  tab$px <- tab$frac$nosZII + 0.5 * tab$frac$nosZIII
  tab$py <- sqrt(3) / 2 * tab$frac$nosZIII
  tab$Total <- rowSums(tab[, need])
  tab$Dominant <- need[max.col(as.matrix(tab[, need]), ties.method = "first")]

  # 用 BLAST 长表补充门信息
  phy <- NULL
  for (b in need) {
    tl <- x$tax_long[[b]]
    if (is.null(tl)) next
    pp <- tax_parse(tl$lineage)
    phy <- rbind(phy, data.frame(Taxon = tax_parse(tl$lineage)[[RANK_LABEL[[level]]]],
                                 Phylum = pp$Phylum, stringsAsFactors = FALSE))
  }
  if (!is.null(phy)) {
    phy <- phy[!duplicated(phy$Taxon), , drop = FALSE]
    tab <- merge(tab, phy, by = "Taxon", all.x = TRUE)
  } else {
    tab$Phylum <- "unassigned"
  }
  tab$Phylum[is.na(tab$Phylum)] <- "unassigned"

  tri <- data.frame(
    x = c(0, 1, 0.5, 0), y = c(0, 0, sqrt(3) / 2, 0)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = tri, ggplot2::aes(x, y),
                          fill = NA, colour = "black", linewidth = 0.4) +
    ggplot2::geom_segment(data = data.frame(
      x = c(0, 0.5, 1), y = c(0, sqrt(3) / 2, 0),
      xend = c(0.5, 0.25, 0.75), yend = c(sqrt(3) / 2, 0, sqrt(3) / 4)
    ), ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey80", linewidth = 0.2)
  if (colour_by == "phylum") {
    p <- p + ggplot2::geom_point(data = tab,
                                 ggplot2::aes(px, py, fill = Phylum, size = Total),
                                 shape = 21, colour = "grey20", stroke = 0.25, alpha = 0.9) +
      ggplot2::scale_fill_manual(values = nosz_pal(length(unique(tab$Phylum)), palette))
  } else {
    p <- p + ggplot2::geom_point(data = tab,
                                 ggplot2::aes(px, py, fill = log10(Total), size = Total),
                                 shape = 21, colour = "grey20", stroke = 0.25, alpha = 0.9) +
      ggplot2::scale_fill_viridis_c(name = expression(log[10] * " total"))
  }
  if (label) {
    p <- p + ggrepel::geom_text_repel(data = tab, ggplot2::aes(px, py, label = Taxon),
                                      size = 1.7, max.overlaps = 20, seed = 1,
                                      min.segment.length = 0.1,
                                      segment.size = 0.15)
  }
  p +
    ggplot2::scale_size_area(max_size = 3.5, name = "Total\nabundance") +
    ggplot2::annotate("text", x = -0.07, y = -0.04, label = "nosZ-I", size = 2.2, hjust = 1) +
    ggplot2::annotate("text", x = 1.07, y = -0.04, label = "nosZ-II", size = 2.2, hjust = 0) +
    ggplot2::annotate("text", x = 0.5, y = sqrt(3) / 2 + 0.06, label = "nosZ-III", size = 2.2) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nosz(border = FALSE, legend_position = "right") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
}

#' Per-taxon bar charts by group
#'
#' A small-multiple view: one panel per taxon, with group means and error bars,
#' individual samples shown as points.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param top_n Number of taxa.
#' @param value `"rpkm"` or `"relative"`.
#' @param palette Palette name.
#' @param ncol Number of columns in the panel grid.
#' @return A ggplot object.
#' @export
plot_top_taxa <- function(x, branch, level = "genus", top_n = 12,
                          value = c("rpkm", "relative"), palette = "nosz",
                          ncol = 4) {
  stopifnot(inherits(x, "noszref"))
  value <- match.arg(value)
  m <- nosz_tax_matrix(x, branch, level)
  if (nrow(m) > top_n) {
    m <- m[order(rowMeans(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  }
  if (value == "relative") m <- rel_abund(m) * 100
  long <- matrix_to_long(m, x$design)
  long$Taxon <- factor(long$Taxon, levels = rownames(m))
  long$SampleID <- factor(long$SampleID, levels = intersect(x$design$SampleID, colnames(m)))
  agg <- long %>% dplyr::group_by(Taxon, Group) %>%
    dplyr::summarise(mean = mean(Abundance, na.rm = TRUE),
                     sd = stats::sd(Abundance, na.rm = TRUE), .groups = "drop")
  ggplot2::ggplot(agg, ggplot2::aes(Group, mean, fill = Group)) +
    ggplot2::geom_col(width = 0.62, colour = "black", linewidth = 0.22) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
                           width = 0.16, linewidth = 0.22) +
    ggplot2::geom_point(data = long, ggplot2::aes(Group, Abundance),
                        inherit.aes = FALSE, shape = 21, fill = "white",
                        size = 0.9, stroke = 0.2,
                        position = ggplot2::position_jitter(width = 0.08)) +
    ggplot2::facet_wrap(~Taxon, scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = group_palette(long$Group, palette)) +
    ggplot2::labs(x = NULL,
                  y = if (value == "relative") "Relative abundance (%)" else "RPKM",
                  title = paste0(branch, " - top ", nrow(m), " ", RANK_LABEL[[level]], " by group")) +
    theme_nosz(x_angle = 30, legend_position = "none")
}

#' Correlation heatmap between abundant taxa
#'
#' Spearman correlation of the most abundant genera (or any rank), ordered by
#' hierarchical clustering and drawn as a triangular tile plot.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param top_n Number of taxa.
#' @param value `"rpkm"` or `"relative"`.
#' @param method Correlation method.
#' @param only_sig Logical, blank out non-significant tiles.
#' @param palette Diverging palette name.
#' @return A ggplot object.
#' @export
plot_tax_corr <- function(x, branch, level = "genus", top_n = 20,
                          value = c("relative", "rpkm"), method = "spearman",
                          only_sig = FALSE, palette = "RdYlBu") {
  stopifnot(inherits(x, "noszref"))
  value <- match.arg(value)
  m <- nosz_tax_matrix(x, branch, level)
  if (nrow(m) > top_n) {
    m <- m[order(rowMeans(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  }
  if (value == "relative") m <- rel_abund(m)
  # 去掉全零行，避免相关系数 NaN
  m <- m[rowSums(m) > 0, , drop = FALSE]
  if (nrow(m) < 3) nosz_abort("可用分类单元少于 3 个，无法计算相关矩阵。")
  ct <- suppressWarnings(stats::cor(t(m), method = method))
  ct[is.na(ct)] <- 0
  pv <- matrix(1, nrow(ct), ncol(ct), dimnames = dimnames(ct))
  for (i in seq_len(nrow(ct))) for (j in seq_len(ncol(ct))) {
    if (i < j) {
      pv[i, j] <- suppressWarnings(stats::cor.test(m[i, ], m[j, ],
                                                   method = method)$p.value)
    }
  }
  ord <- try(stats::hclust(stats::as.dist(1 - ct))$order, silent = TRUE)
  if (!inherits(ord, "try-error")) {
    ct <- ct[ord, ord]; pv <- pv[ord, ord]
  }
  df <- data.frame(
    Taxa1 = rep(rownames(ct), times = ncol(ct)),
    Taxa2 = rep(colnames(ct), each = nrow(ct)),
    r = as.numeric(ct), p = as.numeric(pv)
  )
  df$Taxa1 <- factor(df$Taxa1, levels = rownames(ct))
  df$Taxa2 <- factor(df$Taxa2, levels = rev(rownames(ct)))
  df$lab <- ifelse(df$p < 0.001, "***", ifelse(df$p < 0.01, "**",
                                               ifelse(df$p < 0.05, "*", "")))
  if (only_sig) df$r[df$p >= 0.05] <- NA
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, palette))(101)
  ggplot2::ggplot(df, ggplot2::aes(Taxa1, Taxa2, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(r), "", lab)),
                       size = 1.8, vjust = 0.75) +
    ggplot2::scale_fill_gradientn(colours = pal, limits = c(-1, 1),
                                  na.value = "grey92", name = paste0("r (", method, ")")) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0(branch, " - ", RANK_LABEL[[level]], " co-occurrence (",
                                 method, ")")) +
    theme_nosz(x_angle = 90, legend_position = "right", grid = FALSE) +
    ggplot2::coord_equal()
}
