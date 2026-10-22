# ---------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------

#' Differential abundance between two groups
#'
#' Computes per-taxon group means, log2 fold change, Welch/Wilcoxon p values
#' (BH adjusted) and, when at least three replicates are available, a 95%
#' confidence interval of the mean difference. Optionally a DESeq2 negative
#' binomial test is run on the count table.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param group1,group2 Groups to compare; default the first two groups of the
#'   design.
#' @param method `"wilcox"`, `"t.test"` or `"deseq2"`.
#' @param value `"rpkm"` (default) or `"relative"`.
#' @param min_prevalence Minimum fraction of samples with non-zero abundance.
#' @param pseudo Pseudocount added before log transformation.
#' @return A data frame with columns `Taxon`, `mean1`, `mean2`, `log2FC`,
#'   `diff`, `ci_low`, `ci_high`, `p`, `padj`, `direction`, `prevalence`.
#' @export
calc_diff_taxa <- function(x, branch, level = "genus", group1 = NULL,
                           group2 = NULL, method = c("wilcox", "t.test", "deseq2"),
                           value = c("rpkm", "relative"), min_prevalence = 0.2,
                           pseudo = 1e-6) {
  stopifnot(inherits(x, "noszref"))
  method <- match.arg(method); value <- match.arg(value)
  m <- nosz_tax_matrix(x, branch, level, value = if (value == "relative") "rpkm" else value)
  if (value == "relative") m <- rel_abund(m) * 100
  d <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
  grp <- unique(d$Group)
  if (length(grp) < 2) nosz_abort("")
  g1 <- group1 %||% grp[1]
  g2 <- group2 %||% grp[2]
  if (!all(c(g1, g2) %in% d$Group)) {
    nosz_abort("", paste(grp, collapse = ", "))
  }
  s1 <- d$SampleID[d$Group == g1]; s2 <- d$SampleID[d$Group == g2]
  m <- m[, c(s1, s2), drop = FALSE]
  m <- m[rowSums(m) > 0, , drop = FALSE]
  if (nrow(m) == 0) nosz_abort("")

  prev <- rowMeans(m > 0)
  keep <- prev >= min_prevalence
  if (sum(keep) < 3) keep <- rep(TRUE, nrow(m))
  m <- m[keep, , drop = FALSE]; prev <- prev[keep]

  m1 <- rowMeans(m[, s1, drop = FALSE]); m2 <- rowMeans(m[, s2, drop = FALSE])
  out <- data.frame(
    Taxon = rownames(m), mean1 = m1, mean2 = m2,
    log2FC = log2((m1 + pseudo) / (m2 + pseudo)),
    diff = m1 - m2, ci_low = NA_real_, ci_high = NA_real_,
    p = NA_real_, prevalence = prev, stringsAsFactors = FALSE
  )

  if (method == "deseq2") {
    cnt <- nosz_tax_matrix(x, branch, level, value = "count")
    cnt <- cnt[rownames(m), c(s1, s2), drop = FALSE]
    cnt <- round(cnt)
    cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]
    ok <- try({
      cd <- data.frame(Group = factor(c(rep(g1, length(s1)), rep(g2, length(s2))),
                                      levels = c(g1, g2)))
      dds <- DESeq2::DESeqDataSetFromMatrix(cnt, cd, ~Group)
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      rr <- as.data.frame(DESeq2::results(dds))
      data.frame(Taxon = rownames(rr), log2FC_ds = rr$log2FoldChange,
                 p_ds = rr$pvalue, padj_ds = rr$padj)
    }, silent = TRUE)
    if (!inherits(ok, "try-error")) {
      out <- merge(out, ok, by = "Taxon", all.x = TRUE)
      out$log2FC <- ifelse(is.na(out$log2FC_ds), out$log2FC, out$log2FC_ds)
      out$p <- out$p_ds
      out$padj <- out$padj_ds
      out$log2FC_ds <- NULL; out$p_ds <- NULL; out$padj_ds <- NULL
    } else {
      method <- "wilcox"
    }
  }

  if (!identical(method, "deseq2") || all(is.na(out$p))) {
    for (i in seq_len(nrow(out))) {
      v1 <- m[i, s1]; v2 <- m[i, s2]
      if (all(c(v1, v2) == 0)) next
      if (method == "t.test") {
        tt <- try(stats::t.test(v1, v2), silent = TRUE)
        if (!inherits(tt, "try-error")) {
          out$p[i] <- tt$p.value
          out$diff[i] <- tt$estimate[1] - tt$estimate[2]
          out$ci_low[i] <- tt$conf.int[1]; out$ci_high[i] <- tt$conf.int[2]
        }
      } else {
        wt <- try(stats::wilcox.test(v1, v2, exact = FALSE), silent = TRUE)
        if (!inherits(wt, "try-error")) out$p[i] <- wt$p.value
        if (length(v1) >= 3 && length(v2) >= 3) {
          tt <- try(stats::t.test(v1, v2), silent = TRUE)
          if (!inherits(tt, "try-error")) {
            out$ci_low[i] <- tt$conf.int[1]; out$ci_high[i] <- tt$conf.int[2]
          }
        }
      }
    }
  }
  out$padj <- stats::p.adjust(out$p, method = "BH")
  out$direction <- ifelse(is.na(out$padj), "ns",
                          ifelse(out$padj < 0.05 & out$log2FC > 0, "Up in first",
                                 ifelse(out$padj < 0.05 & out$log2FC < 0, "Down in first", "ns")))
  out$Group1 <- g1; out$Group2 <- g2
  out$Branch <- branch; out$Level <- level; out$Method <- method
  out <- out[order(out$padj, -abs(out$log2FC)), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "groups") <- c(g1, g2)
  out
}

# 
diff_get <- function(x, res, branch, level, method, ...) {
  if (!is.null(res)) return(res)
  calc_diff_taxa(x, branch = branch, level = level, method = method, ...)
}

diff_cols <- c("Up in first" = "#D7191C", "Down in first" = "#2C7AB6", "ns" = "grey75")

#' Volcano plot of differential taxa
#'
#' @param x A `noszref` object.
#' @param res Optional result of [calc_diff_taxa()].
#' @param branch Branch name (used when `res` is `NULL`).
#' @param level Taxonomic level.
#' @param method Differential test.
#' @param padj_cut,fc_cut Significance thresholds.
#' @param label_n Number of taxa to label.
#' @param palette Colour palette for up/down/ns.
#' @return A ggplot object.
#' @export
plot_diff_volcano <- function(x, res = NULL, branch = "nosZI", level = "genus",
                              method = "wilcox", padj_cut = 0.05, fc_cut = 1,
                              label_n = 10, palette = NULL) {
  stopifnot(inherits(x, "noszref"))
  res <- diff_get(x, res, branch, level, method)
  res$sig <- ifelse(is.na(res$padj), "ns",
                    ifelse(res$padj < padj_cut & abs(res$log2FC) >= fc_cut,
                           ifelse(res$log2FC > 0, "Up in first", "Down in first"), "ns"))
  res$y <- -log10(pmax(res$padj, 1e-300))
  cols <- if (is.null(palette)) diff_cols else
    stats::setNames(nosz_pal(3, palette), names(diff_cols))
  lab <- res[res$sig != "ns", , drop = FALSE]
  lab <- head(lab[order(lab$padj), , drop = FALSE], label_n)
  g <- attr(res, "groups") %||% c(res$Group1[1], res$Group2[1])
  xlab_txt <- bquote(log[2] * " fold change (" * .(g[1]) * " / " * .(g[2]) * ")")
  ggplot2::ggplot(res, ggplot2::aes(log2FC, y, colour = sig)) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = "dashed",
                        linewidth = 0.25, colour = "grey50") +
    ggplot2::geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed",
                        linewidth = 0.25, colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(size = prevalence), alpha = 0.85) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = Taxon),
                             size = 1.8, max.overlaps = 30, seed = 2,
                             min.segment.length = 0.1, segment.size = 0.15,
                             show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_size_continuous(range = c(0.3, 2.4), name = "Prevalence") +
    ggplot2::labs(x = xlab_txt,
                  y = expression(-log[10] * "(adjusted " * italic(P) * ")"),
                  title = paste0(res$Branch[1], " - ", RANK_LABEL[[res$Level[1]]],
                                 " differential abundance")) +
    theme_nosz(grid = TRUE, legend_position = "right")
}

#' MA plot of differential taxa
#'
#' @inheritParams plot_diff_volcano
#' @return A ggplot object.
#' @export
#' @param padj_cut P-value cutoff
plot_diff_ma <- function(x, res = NULL, branch = "nosZI", level = "genus",
                         method = "wilcox", padj_cut = 0.05, label_n = 8) {
  stopifnot(inherits(x, "noszref"))
  res <- diff_get(x, res, branch, level, method)
  res$A <- log10((res$mean1 + res$mean2) / 2 + 1)
  res$M <- res$log2FC
  res$sig <- ifelse(!is.na(res$padj) & res$padj < padj_cut,
                    ifelse(res$M > 0, "Up in first", "Down in first"), "ns")
  lab <- head(res[res$sig != "ns", , drop = FALSE][order(res$padj[res$sig != "ns"]), ],
              label_n)
  ggplot2::ggplot(res, ggplot2::aes(A, M, colour = sig)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey40") +
    ggplot2::geom_point(size = 1.1, alpha = 0.85) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = Taxon), size = 1.8,
                             max.overlaps = 20, seed = 3, min.segment.length = 0.1,
                             segment.size = 0.15, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = diff_cols, name = NULL) +
    ggplot2::labs(x = expression(log[10] * "(mean abundance + 1)"),
                  y = expression(log[2] * " fold change"),
                  title = paste0(res$Branch[1], " - MA plot")) +
    theme_nosz(grid = TRUE, legend_position = "right")
}

#' Forest plot of effect sizes
#'
#' @inheritParams plot_diff_volcano
#' @param top_n Number of taxa shown (ranked by absolute log2FC among
#'   significant taxa, falling back to all taxa).
#' @param show_ci Logical, draw 95% confidence intervals when available.
#' @return A ggplot object.
#' @export
#' @param padj_cut P-value cutoff
plot_diff_forest <- function(x, res = NULL, branch = "nosZI", level = "genus",
                             method = "t.test", padj_cut = 0.05, top_n = 15,
                             show_ci = TRUE) {
  stopifnot(inherits(x, "noszref"))
  res <- diff_get(x, res, branch, level, method)
  sig <- res[!is.na(res$padj) & res$padj < padj_cut, , drop = FALSE]
  if (nrow(sig) < 2) sig <- res
  sig <- head(sig[order(-abs(sig$log2FC)), , drop = FALSE], top_n)
  sig$Taxon <- factor(sig$Taxon, levels = sig$Taxon[order(sig$log2FC)])
  sig$dir <- ifelse(sig$log2FC > 0, "Up in first", "Down in first")
  g <- attr(res, "groups") %||% c(res$Group1[1], res$Group2[1])
  xlab_txt <- bquote(log[2] * " FC (" * .(g[1]) * " / " * .(g[2]) * ")")
  p <- ggplot2::ggplot(sig, ggplot2::aes(log2FC, Taxon, colour = dir)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey40") +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = diff_cols, name = NULL) +
    ggplot2::labs(x = xlab_txt, y = NULL,
                  title = paste0(res$Branch[1], " - effect size (",
                                 RANK_LABEL[[res$Level[1]]], ")")) +
    theme_nosz(grid = TRUE, legend_position = "top")
  if (show_ci && any(!is.na(sig$ci_low))) {
    p <- p + ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high),
                                     height = 0.25, linewidth = 0.3, na.rm = TRUE)
  }
  p
}

#' STAMP-style extended error bar plot
#'
#' Two-sided bar chart of group means with 95% confidence intervals, plus a
#' central panel showing the mean difference for each taxon -- the layout used
#' by STAMP and frequently requested by reviewers.
#'
#' @inheritParams plot_diff_volcano
#' @param top_n Number of taxa to display.
#' @param sort_by `"difference"` or `"p"`.
#' @param palette Palette name for the two groups.
#' @param width_ratio Relative widths of the left, centre and right panels.
#' @return A patchwork object.
#' @export
#' @param padj_cut P-value cutoff
plot_diff_stamp <- function(x, res = NULL, branch = "nosZI", level = "genus",
                            method = "t.test", padj_cut = 0.05, top_n = 12,
                            sort_by = c("difference", "p"), palette = "nosz",
                            width_ratio = c(3, 1.2, 3)) {
  stopifnot(inherits(x, "noszref"))
  sort_by <- match.arg(sort_by)
  res <- diff_get(x, res, branch, level, method, value = "relative")
  res$level1 <- res$mean1; res$level2 <- res$mean2
  sel <- res[!is.na(res$padj) & res$padj < padj_cut, , drop = FALSE]
  if (nrow(sel) < 2) sel <- res
  if (sort_by == "difference") {
    sel <- head(sel[order(-abs(sel$diff)), , drop = FALSE], top_n)
  } else {
    sel <- head(sel[order(sel$padj), , drop = FALSE], top_n)
  }
  g <- attr(res, "groups") %||% c(res$Group1[1], res$Group2[1])
  ord <- sel$Taxon[order(sel$diff)]
  sel$Taxon <- factor(sel$Taxon, levels = ord)
  cols <- group_palette(c(g[1], g[2]), palette)

  #  95% t  STAMP 
  m <- rel_abund(nosz_tax_matrix(x, branch, level)) * 100
  s1 <- intersect(x$design$SampleID[x$design$Group == g[1]], colnames(m))
  s2 <- intersect(x$design$SampleID[x$design$Group == g[2]], colnames(m))
  ci_mean <- function(v) {
    v <- v[!is.na(v)]; n <- length(v)
    if (n < 2) return(c(mean(v), mean(v)))
    se <- stats::sd(v) / sqrt(n)
    mean(v) + c(-1, 1) * stats::qt(0.975, n - 1) * se
  }
  c1 <- t(vapply(as.character(sel$Taxon), function(t) ci_mean(m[t, s1]), numeric(2)))
  c2 <- t(vapply(as.character(sel$Taxon), function(t) ci_mean(m[t, s2]), numeric(2)))
  long <- rbind(
    data.frame(Taxon = sel$Taxon, Group = g[1], Mean = sel$mean1,
               Lo = c1[, 1], Hi = c1[, 2]),
    data.frame(Taxon = sel$Taxon, Group = g[2], Mean = sel$mean2,
               Lo = c2[, 1], Hi = c2[, 2])
  )
  long$Group <- factor(long$Group, levels = c(g[1], g[2]))

  p_left <- ggplot2::ggplot(long, ggplot2::aes(Mean, Taxon, fill = Group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72),
                      width = 0.66, colour = "black", linewidth = 0.2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = pmax(Lo, 0), xmax = Hi),
                            height = 0.18, linewidth = 0.25, na.rm = TRUE,
                            position = ggplot2::position_dodge(width = 0.72)) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "Mean relative abundance (%)", y = NULL) +
    theme_nosz(legend_position = "top", grid = TRUE)

  p_mid <- ggplot2::ggplot(sel, ggplot2::aes(diff, Taxon)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey40") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = diff, yend = Taxon),
                          linewidth = 0.3, colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(colour = diff > 0), size = 1.3, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "#2C7AB6")) +
    ggplot2::labs(x = "Difference (%)", y = NULL) +
    theme_nosz(grid = FALSE, legend_position = "none") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())

  p_right <- ggplot2::ggplot(sel, ggplot2::aes(-log10(pmax(padj, 1e-300)), Taxon)) +
    ggplot2::geom_col(width = 0.6, fill = "grey55", colour = "black", linewidth = 0.15) +
    ggplot2::geom_vline(xintercept = -log10(padj_cut), linetype = "dashed",
                        linewidth = 0.25, colour = "#D7191C") +
    ggplot2::labs(x = expression(-log[10] * "(" * italic(P)[adj] * ")"), y = NULL) +
    theme_nosz(grid = TRUE, legend_position = "none") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())

  patchwork::wrap_plots(p_left, p_mid, p_right, nrow = 1,
                        widths = width_ratio) +
    patchwork::plot_annotation(
      title = paste0(res$Branch[1], " - STAMP-style comparison (",
                     RANK_LABEL[[res$Level[1]]], ")"),
      subtitle = paste0(g[1], " vs ", g[2], "; ", nrow(sel), " taxa shown"))
}

#' Heatmap of differential taxa
#'
#' @inheritParams plot_diff_volcano
#' @param top_n Number of differential taxa.
#' @param scale `"row"` or `"none"`.
#' @param show_stars Logical, annotate tiles with significance stars.
#' @return A `pheatmap` object.
#' @export
#' @param padj_cut P-value cutoff
plot_diff_heatmap <- function(x, res = NULL, branch = "nosZI", level = "genus",
                              method = "wilcox", padj_cut = 0.05, top_n = 20,
                              scale = "row", show_stars = TRUE) {
  stopifnot(inherits(x, "noszref"))
  res <- diff_get(x, res, branch, level, method)
  sel <- res[!is.na(res$padj) & res$padj < padj_cut, , drop = FALSE]
  if (nrow(sel) < 2) sel <- res
  sel <- head(sel[order(sel$padj), , drop = FALSE], top_n)
  m <- nosz_tax_matrix(x, branch, level, value = "rpkm")
  m <- m[intersect(sel$Taxon, rownames(m)), , drop = FALSE]
  if (nrow(m) == 0) nosz_abort(" RPKM ")
  m <- log10(m + 1)
  m <- m[, order_samples(colnames(m), x$design), drop = FALSE]
  ann <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
  rownames(ann) <- ann$SampleID
  ann$SampleID <- NULL; ann[] <- lapply(ann, as.character)
  stars <- NULL
  if (show_stars) {
    st <- p_star(sel$padj[match(rownames(m), sel$Taxon)])
    stars <- matrix(st, nrow = nrow(m), ncol = ncol(m))
    rownames(stars) <- rownames(m); colnames(stars) <- colnames(m)
  }
  pheatmap::pheatmap(
    m, scale = if (scale == "none") "none" else scale,
    cluster_rows = nrow(m) > 2, cluster_cols = FALSE,
    annotation_col = ann, border_colour = "white", border_width = 0.6,
    display_numbers = if (show_stars) stars else FALSE,
    fontsize = 7, fontsize_number = 6, angle_col = 45,
    color = grDevices::colorRampPalette(c("#2C7AB6", "white", "#D7191C"))(101),
    main = paste0(res$Branch[1], " - top differential ",
                  RANK_LABEL[[res$Level[1]]], " (", res$Group1[1], " vs ",
                  res$Group2[1], ")")
  )
}
