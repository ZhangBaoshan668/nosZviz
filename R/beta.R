# ---------------------------------------------------------------------------
# beta 
# ---------------------------------------------------------------------------

#  design 
beta_matrix <- function(x, branch, distance = "bray_curtis") {
  mats <- nosz_need(x, "beta", branch)
  if (is.null(mats[[distance]])) {
    nosz_abort(" ", branch, "  `", distance, "`",
               paste(names(mats), collapse = ", "))
  }
  m <- mats[[distance]]
  m <- (m + t(m)) / 2
  diag(m) <- 0
  smp <- intersect(x$design$SampleID, rownames(m))
  if (length(smp) < 2) smp <- rownames(m)
  m[smp, smp, drop = FALSE]
}

# 
beta_design <- function(x, m) {
  d <- x$design[x$design$SampleID %in% rownames(m), , drop = FALSE]
  d <- d[match(rownames(m), d$SampleID), , drop = FALSE]
  d$SampleID <- factor(d$SampleID, levels = rownames(m))
  d$Group <- factor(d$Group, levels = unique(d$Group))
  d
}

#' PERMANOVA, ANOSIM, MRPP and beta-dispersion tests
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param distance Distance metric, one of those available in the object
#'   (`"bray_curtis"`, `"jaccard"`, `"euclidean"`, `"manhatten"`, and the
#'   `"_binary"` variants).
#' @param permutations Number of permutations.
#' @return A data frame with method, statistic and p value.
#' @export
calc_permanova <- function(x, branch, distance = "bray_curtis",
                           permutations = 999) {
  stopifnot(inherits(x, "noszref"))
  m <- beta_matrix(x, branch, distance)
  d <- beta_design(x, m)
  if (nlevels(d$Group) < 2) {
    return(data.frame(Method = "PERMANOVA", Statistic = NA_real_, p = NA_real_,
                      Note = "only one group", stringsAsFactors = FALSE))
  }
  dd <- stats::as.dist(m)
  out <- list()
  ad <- try(vegan::adonis2(dd ~ Group, data = d, permutations = permutations),
            silent = TRUE)
  if (!inherits(ad, "try-error")) {
    tb <- as.data.frame(ad)
    out[[1]] <- data.frame(Method = "PERMANOVA",
                           Statistic = tb$R2[1],
                           p = tb$`Pr(>F)`[1],
                           Note = paste0("R2 = ", formatC(tb$R2[1], 3, format = "f")),
                           stringsAsFactors = FALSE)
  }
  an <- try(vegan::anosim(dd, d$Group, permutations = permutations), silent = TRUE)
  if (!inherits(an, "try-error")) {
    out[[length(out) + 1]] <- data.frame(
      Method = "ANOSIM", Statistic = unname(an$statistic),
      p = an$signif, Note = paste0("R = ", formatC(an$statistic, 3, format = "f")),
      stringsAsFactors = FALSE)
  }
  mp <- try(vegan::mrpp(dd, d$Group, permutations = permutations), silent = TRUE)
  if (!inherits(mp, "try-error")) {
    out[[length(out) + 1]] <- data.frame(
      Method = "MRPP", Statistic = unname(mp$A), p = mp$Pvalue,
      Note = paste0("A = ", formatC(mp$A, 3, format = "f")),
      stringsAsFactors = FALSE)
  }
  bd <- try(vegan::betadisper(dd, d$Group), silent = TRUE)
  if (!inherits(bd, "try-error")) {
    a <- try(stats::anova(bd), silent = TRUE)
    if (!inherits(a, "try-error")) {
      out[[length(out) + 1]] <- data.frame(
        Method = "Beta-dispersion", Statistic = a$`F value`[1],
        p = a$`Pr(>F)`[1], Note = paste0("F = ", formatC(a$`F value`[1], 3, format = "f")),
        stringsAsFactors = FALSE)
    }
  }
  if (length(out) == 0) nosz_abort("")
  res <- do.call(rbind, out)
  res$Distance <- distance
  res
}

# PCoA / NMDS / PCA
.ordination <- function(x, branch, method = c("pcoa", "nmds", "pca"),
                        distance = "bray_curtis") {
  method <- match.arg(method)
  if (method == "pca") {
    otu <- nosz_otu_matrix(x, branch, rarefied = FALSE)
    if (is.null(otu)) otu <- nosz_otu_matrix(x, branch, rarefied = TRUE)
    otu <- otu[, colSums(otu) > 0, drop = FALSE]
    hel <- vegan::decostand(t(otu + 0.5), method = "hellinger")
    pc <- stats::prcomp(hel, center = TRUE, scale. = FALSE)
    sc <- as.data.frame(pc$x[, 1:2, drop = FALSE])
    ve <- (pc$sdev^2 / sum(pc$sdev^2))[1:2] * 100
    colnames(sc) <- c("Axis1", "Axis2")
    return(list(scores = sc, var = ve, method = "PCA (Hellinger)",
                xlab = paste0("PC1 (", formatC(ve[1], 1, format = "f"), "%)"),
                ylab = paste0("PC2 (", formatC(ve[2], 1, format = "f"), "%)"),
                note = NULL))
  }
  m <- beta_matrix(x, branch, distance)
  d <- stats::as.dist(m)
  if (method == "pcoa") {
    fit <- stats::cmdscale(d, k = 2, eig = TRUE)
    pos <- fit$eig[fit$eig > 0]
    ve <- fit$eig[1:2] / sum(pos) * 100
    sc <- as.data.frame(fit$points[, 1:2, drop = FALSE])
    colnames(sc) <- c("Axis1", "Axis2")
    return(list(scores = sc, var = ve, method = "PCoA",
                xlab = paste0("PCoA1 (", formatC(ve[1], 1, format = "f"), "%)"),
                ylab = paste0("PCoA2 (", formatC(ve[2], 1, format = "f"), "%)"),
                note = NULL))
  }
  fit <- try(vegan::metaMDS(d, k = 2, trymax = 200, autotransform = FALSE,
                            trace = FALSE, plot = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) nosz_abort("NMDS  PCoA")
  sc <- as.data.frame(fit$points[, 1:2, drop = FALSE])
  colnames(sc) <- c("Axis1", "Axis2")
  list(scores = sc, var = NA, method = "NMDS",
       xlab = "NMDS1", ylab = "NMDS2",
       note = paste0("Stress = ", formatC(fit$stress, 3, format = "f")))
}

#' Beta diversity ordination (PCoA / NMDS / PCA)
#'
#' Draws samples in a two-dimensional ordination with group confidence
#' ellipses, optional sample labels and a PERMANOVA annotation.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param method `"pcoa"` (default), `"nmds"` or `"pca"`.
#' @param distance Distance metric (ignored for `"pca"`).
#' @param ellipse Logical, draw 95% confidence ellipses (needs >= 3 samples
#'   per group).
#' @param label Logical, label samples with [ggrepel].
#' @param permanova Logical, run PERMANOVA and annotate the subtitle.
#' @param permutations Permutations for PERMANOVA.
#' @param palette Palette name.
#' @param hull Logical, draw convex hulls instead of ellipses.
#' @param star Logical, connect samples to group centroid.
#' @return A ggplot object.
#' @export
plot_ordination <- function(x, branch, method = c("pcoa", "nmds", "pca"),
                            distance = "bray_curtis", ellipse = TRUE,
                            label = TRUE, permanova = TRUE,
                            permutations = 999, palette = "nosz",
                            hull = FALSE, star = FALSE) {
  stopifnot(inherits(x, "noszref"))
  method <- match.arg(method)
  ord <- .ordination(x, branch, method, distance)
  sc <- ord$scores
  sc$SampleID <- rownames(sc)
  d <- x$design[x$design$SampleID %in% sc$SampleID, , drop = FALSE]
  sc <- merge(sc, d, by = "SampleID", all.x = TRUE)
  sc$SampleID <- factor(sc$SampleID, levels = d$SampleID)
  sc$Group <- factor(sc$Group, levels = unique(d$Group))
  cols <- group_palette(sc$Group, palette)

  sub <- ord$note
  if (permanova && method != "pca" && nlevels(droplevels(sc$Group)) >= 2) {
    pv <- try(calc_permanova(x, branch, distance, permutations), silent = TRUE)
    if (!inherits(pv, "try-error")) {
      r <- pv[pv$Method == "PERMANOVA", , drop = FALSE]
      if (nrow(r) > 0 && !is.na(r$p[1])) {
        sub <- paste0(sub %||% "", ifelse(is.null(sub), "", "; "),
                      "PERMANOVA R2 = ", formatC(r$Statistic[1], 3, format = "f"),
                      ", p = ", fmt_p(r$p[1]))
      }
    }
  }

  p <- ggplot2::ggplot(sc, ggplot2::aes(Axis1, Axis2, colour = Group, fill = Group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2,
                        colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2,
                        colour = "grey70") +
    ggplot2::geom_point(size = 1.9, shape = 21, stroke = 0.3, alpha = 0.95) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = ord$xlab, y = ord$ylab,
                  subtitle = if (nzchar(sub %||% "")) sub else NULL,
                  title = paste0(branch, " - ", ord$method,
                                 if (method != "pca") paste0(" (", distance, ")") else "")) +
    theme_nosz(grid = FALSE, legend_position = "right")
  if (ellipse && all(table(sc$Group) >= 3)) {
    p <- p + ggplot2::stat_ellipse(ggplot2::aes(colour = Group), level = 0.95,
                                   linewidth = 0.35, show.legend = FALSE,
                                   type = "t")
  }
  if (hull) {
    hull_df <- sc %>% dplyr::group_by(Group) %>%
      dplyr::slice(grDevices::chull(Axis1, Axis2)) %>% dplyr::ungroup()
    p <- p + ggplot2::geom_polygon(data = hull_df, ggplot2::aes(fill = Group),
                                   alpha = 0.1, linetype = "dashed",
                                   linewidth = 0.25, show.legend = FALSE)
  }
  if (star) {
    cen <- sc %>% dplyr::group_by(Group) %>%
      dplyr::summarise(Axis1 = mean(Axis1), Axis2 = mean(Axis2), .groups = "drop")
    p <- p + ggplot2::geom_segment(data = merge(sc, cen, by = "Group",
                                                suffixes = c("", ".c")),
                                   ggplot2::aes(xend = Axis1.c, yend = Axis2.c),
                                   linewidth = 0.2, alpha = 0.55) +
      ggplot2::geom_point(data = cen, ggplot2::aes(Axis1, Axis2, fill = Group),
                          shape = 23, size = 2.2, stroke = 0.3)
  }
  if (label) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = SampleID), size = 1.9,
                                      seed = 1, min.segment.length = 0.15,
                                      segment.size = 0.15, show.legend = FALSE)
  }
  p
}

#' Heatmap of a beta diversity distance matrix
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param distance Distance metric.
#' @param annotation Logical, add group annotation bar.
#' @param cluster Logical, cluster samples.
#' @param palette Sequential palette name.
#' @param numbers Logical, show distance values in cells.
#' @return A `pheatmap` object.
#' @export
plot_beta_heatmap <- function(x, branch, distance = "bray_curtis",
                              annotation = TRUE, cluster = TRUE,
                              palette = "YlOrRd", numbers = TRUE) {
  stopifnot(inherits(x, "noszref"))
  m <- beta_matrix(x, branch, distance)
  ann <- NULL
  if (annotation) {
    ann <- x$design[x$design$SampleID %in% rownames(m), , drop = FALSE]
    rownames(ann) <- ann$SampleID
    ann$SampleID <- NULL
    ann[] <- lapply(ann, as.character)
  }
  pheatmap::pheatmap(
    m, cluster_rows = cluster, cluster_cols = cluster,
    annotation_col = ann, annotation_row = ann,
    border_colour = "white", border_width = 0.8,
    display_numbers = numbers, number_format = "%.2f",
    fontsize = 7, fontsize_number = 5.5,
    color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, palette))(101),
    angle_col = 45, main = paste0(branch, " - ", distance, " distance")
  )
}

#' UPGMA dendrogram with group annotation
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param distance Distance metric.
#' @param method Linkage method, see [stats::hclust()].
#' @param palette Palette name.
#' @param label_size Label size.
#' @return A ggplot object.
#' @export
plot_beta_dendrogram <- function(x, branch, distance = "bray_curtis",
                                 method = "average", palette = "nosz",
                                 label_size = 2.1) {
  stopifnot(inherits(x, "noszref"))
  m <- beta_matrix(x, branch, distance)
  d <- beta_design(x, m)
  hc <- stats::hclust(stats::as.dist(m), method = method)
  lay <- ggraph::create_layout(hc, layout = "dendrogram")
  #  layout  merge 
  lay$Group <- d$Group[match(lay$label, as.character(d$SampleID))]
  cols <- group_palette(d$Group, palette)
  ggraph::ggraph(lay) +
    ggraph::geom_edge_elbow(linewidth = 0.4, colour = "grey35") +
    ggraph::geom_node_point(ggplot2::aes(filter = leaf, colour = Group), size = 1.6) +
    ggraph::geom_node_text(ggplot2::aes(filter = leaf, label = label),
                           size = label_size, angle = 90, hjust = 1.15,
                           vjust = 0.5) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = paste0("Distance (", method, " linkage)"),
                  title = paste0(branch, " - UPGMA on ", distance)) +
    theme_nosz(border = FALSE, legend_position = "right") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
}

#' Beta dispersion (multivariate homogeneity) boxplot
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param distance Distance metric.
#' @param palette Palette name.
#' @return A ggplot object.
#' @export
plot_beta_dispersion <- function(x, branch, distance = "bray_curtis",
                                 palette = "nosz") {
  stopifnot(inherits(x, "noszref"))
  m <- beta_matrix(x, branch, distance)
  d <- beta_design(x, m)
  if (nlevels(d$Group) < 2) nosz_abort(" beta dispersion")
  bd <- vegan::betadisper(stats::as.dist(m), d$Group)
  df <- data.frame(SampleID = names(bd$distances),
                   Distance = as.numeric(bd$distances),
                   Group = bd$group)
  a <- stats::anova(bd)
  p_lab <- paste0("ANOVA F = ", formatC(a$`F value`[1], 3, format = "f"),
                  ", p = ", fmt_p(a$`Pr(>F)`[1]))
  ggplot2::ggplot(df, ggplot2::aes(Group, Distance, fill = Group)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA,
                          colour = "black", linewidth = 0.3, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.1, size = 1.1, shape = 21, fill = "white",
                         colour = "black", stroke = 0.25) +
    ggplot2::scale_fill_manual(values = group_palette(df$Group, palette)) +
    ggplot2::labs(x = NULL, y = "Distance to group centroid",
                  subtitle = p_lab,
                  title = paste0(branch, " - beta dispersion (", distance, ")")) +
    theme_nosz(x_angle = 30, legend_position = "none")
}

#' Multi-distance ordination panel
#'
#' Runs PCoA on several distance metrics and assembles the panels with a
#' shared legend, a compact way to show that the grouping is robust to the
#' choice of metric.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param distances Distance metrics to compare.
#' @param palette Palette name.
#' @param ncol Columns of the panel grid.
#' @return A patchwork object.
#' @export
plot_beta_panel <- function(x, branch,
                            distances = c("bray_curtis", "jaccard",
                                          "euclidean", "manhatten"),
                            palette = "nosz", ncol = 2) {
  stopifnot(inherits(x, "noszref"))
  avail <- names(nosz_need(x, "beta", branch))
  distances <- intersect(distances, avail)
  if (length(distances) == 0) nosz_abort("")
  plist <- lapply(distances, function(dd) {
    plot_ordination(x, branch, method = "pcoa", distance = dd,
                    permanova = FALSE, palette = palette) +
      ggplot2::labs(title = dd, subtitle = NULL) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(size = 8))
  })
  patchwork::wrap_plots(plist, ncol = ncol) +
    patchwork::plot_annotation(tag_levels = "A",
                               title = paste0(branch, " - PCoA across distance metrics"))
}

#' Sample-sample correlation heatmap
#'
#' Spearman correlation between samples based on OTU relative abundance,
#' shown as a lower triangle with hierarchical ordering.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param method Correlation method.
#' @param rarefied Use the rarefied OTU table.
#' @param triangle Logical, show only the lower triangle.
#' @param palette Diverging palette.
#' @return A ggplot object.
#' @export
plot_sample_corr <- function(x, branch, method = "spearman", rarefied = TRUE,
                             triangle = TRUE, palette = "RdYlBu") {
  stopifnot(inherits(x, "noszref"))
  otu <- try(nosz_otu_matrix(x, branch, rarefied = rarefied), silent = TRUE)
  if (inherits(otu, "try-error")) otu <- nosz_otu_matrix(x, branch, rarefied = FALSE)
  otu <- otu[, colSums(otu) > 0, drop = FALSE]
  if (ncol(otu) < 2) nosz_abort("")
  r <- stats::cor(otu, method = method)
  ord <- try(stats::hclust(stats::as.dist(1 - r))$order, silent = TRUE)
  if (!inherits(ord, "try-error")) r <- r[ord, ord]
  df <- data.frame(
    S1 = rep(rownames(r), times = ncol(r)),
    S2 = rep(colnames(r), each = nrow(r)),
    r = as.numeric(r)
  )
  df$S1 <- factor(df$S1, levels = rownames(r))
  df$S2 <- factor(df$S2, levels = rev(rownames(r)))
  if (triangle) df$r[as.integer(df$S1) < as.integer(df$S2)] <- NA
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, palette))(101)
  ggplot2::ggplot(df, ggplot2::aes(S1, S2, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(r), "", formatC(r, 2, format = "f"))),
                       size = 2) +
    ggplot2::scale_fill_gradientn(colours = pal, na.value = "white",
                                  name = paste0(method, " r")) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0(branch, " - sample correlation (", method, ")")) +
    theme_nosz(x_angle = 45, legend_position = "right", grid = FALSE) +
    ggplot2::coord_equal()
}
