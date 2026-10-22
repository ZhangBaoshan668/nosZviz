#'
# NosZRef 
#  -  RDA / CCA / db-RDA
#  - envfit  PCoA 
#  - LEfSe Kruskal-Wallis + Wilcoxon + 
#  - IndVal + 
#  -  OTU UpSet
#  - alpha   
#  - -
#  -  UPGMA 
#'

# 
prepare_env <- function(x, env) {
  if (is.null(env)) return(NULL)
  if (is.character(env)) {
    env <- read_delim_safe(env, delim = if (grepl(",", readLines(env, 1))) "," else "\t")
  }
  if (!is.data.frame(env)) return(NULL)
  env <- as.data.frame(env, stringsAsFactors = FALSE)
  nm <- tolower(colnames(env))
  id_i <- which(nm %in% c("sampleid", "sample", "sample_id", "id"))[1]
  if (is.na(id_i)) id_i <- 1
  colnames(env)[id_i] <- "SampleID"
  env$SampleID <- as.character(env$SampleID)
  num <- vapply(env, is.numeric, logical(1))
  keep <- colnames(env)[num & colnames(env) != "SampleID"]
  keep <- keep[vapply(env[keep], function(v) stats::sd(v, na.rm = TRUE) > 0,
                      logical(1))]
  if (length(keep) == 0) return(NULL)
  env <- env[, c("SampleID", keep), drop = FALSE]
  env <- env[env$SampleID %in% x$design$SampleID, , drop = FALSE]
  if (nrow(env) < 3) return(NULL)
  env
}

#' Constrained ordination of nosZ taxa with environmental variables
#'
#' nosZ RDA / CCA / db-RDA
#'  vegan  rda()/cca()/capscale()
#'  Hellinger 
#'  P ****
#'
#' @param x A `noszref` object.
#' @param branch Branch name (one branch at a time).
#' @param level Taxonomic level used as response matrix.
#' @param env Environmental data: data frame with a `SampleID` column and
#'   numeric variables, or path to such a table.
#' @param method `"RDA"`Hellinger , `"CCA"` 
#'   `"dbRDA"`.
#' @param top_n Top taxa (by mean abundance) entering the model.
#' @param dist Distance for dbRDA, e.g. `"bray"`.
#' @param permutations Permutations for the model test.
#' @param show_species Show top species scores as text.
#' @param palette Palette for groups.
#' @return A ggplot object; the fitted model is attached as attribute
#'   `"model"`.
#' @export
plot_ordination_constrained <- function(x, branch, level = "genus", env,
                                        method = c("RDA", "CCA", "dbRDA"),
                                        top_n = 30, dist = "bray",
                                        permutations = 999, show_species = TRUE,
                                        palette = "nosz") {
  stopifnot(inherits(x, "noszref"))
  method <- match.arg(method)
  env <- prepare_env(x, env)
  if (is.null(env)) nosz_abort(" SampleID ")
  m <- nosz_tax_matrix(x, branch, level, value = "rpkm")
  smp <- intersect(colnames(m), env$SampleID)
  if (length(smp) < 4) nosz_abort(" 4 ")
  m <- m[, smp, drop = FALSE]
  env <- env[match(smp, env$SampleID), , drop = FALSE]
  rownames(env) <- env$SampleID
  env$SampleID <- NULL
  # top_n   x vegan =
  if (nrow(m) > top_n) m <- m[order(rowMeans(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  ms <- t(m)

  if (method == "RDA") {
    hel <- vegan::decostand(ms, "hellinger")
    ord <- vegan::rda(hel, env)
  } else if (method == "CCA") {
    rel <- vegan::decostand(ms, "total")
    ord <- vegan::cca(rel, env)
  } else {
    ord <- vegan::capscale(vegan::vegdist(ms, method = dist) ~ ., data = env)
  }
  pv <- tryCatch({
    an <- vegan::anova.cca(ord, permutations = permutations)
    sprintf("P = %s", fmt_p(an[1, "Pr(>F)"]))
  }, error = function(e) "P = NA")
  tot <- tryCatch({
    s <- summary(ord)
    sprintf("%.1f%% variance constrained", 100 * s$constr.chi / s$tot.chi)
  }, error = function(e) "variance NA")

  sc <- vegan::scores(ord, display = c("sites", "species", "biplot"),
                      choices = 1:2, scaling = 2)
  si_m <- sc$sites
  #  0  1  0 
  if (is.null(si_m) || ncol(si_m) < 2) {
    ax2 <- rep(0, nrow(si_m))
    si_m <- cbind(si_m, ax2)
    colnames(si_m)[2] <- paste0(colnames(si_m)[1], ".b")
  }
  xcol <- colnames(si_m)[1]
  ycol <- colnames(si_m)[2]
  si <- as.data.frame(si_m)
  si$SampleID <- rownames(si)
  si <- merge(si, x$design[, c("SampleID", "Group")], by = "SampleID", all.x = TRUE)
  cols <- group_palette(si$Group, palette)
  ev <- tryCatch(vegan::eigenvals(ord)[1:2], error = function(e) c(NA, NA))
  xl <- sprintf("RDA%d (%.1f%%)", 1:2, 100 * ev / sum(vegan::eigenvals(ord)))[1]
  yl <- sprintf("RDA%d (%.1f%%)", 1:2, 100 * ev / sum(vegan::eigenvals(ord)))[2]
  if (method == "CCA") { xl <- sub("RDA", "CCA", xl); yl <- sub("RDA", "CCA", yl) }
  if (method == "dbRDA") { xl <- sub("RDA", "CAP", xl); yl <- sub("RDA", "CAP", yl) }

  p <- ggplot2::ggplot(si, ggplot2::aes(.data[[xcol]], .data[[ycol]],
                                        colour = Group)) +
    ggplot2::geom_point(size = 2, alpha = 0.95) +
    ggplot2::stat_ellipse(ggplot2::aes(fill = Group), geom = "polygon",
                          level = 0.68, alpha = 0.12, colour = NA,
                          show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = xl, y = yl,
                  title = paste0(branch, " ", method, " - ", level, " x environment"),
                  subtitle = paste(tot, pv, sep = "; ")) +
    theme_nosz()

  #  0.8vegan  biplot
  bp <- as.data.frame(sc$biplot)
  if (!is.null(bp) && nrow(bp) > 0 && ncol(bp) >= 2) {
    r <- max(abs(c(si[[xcol]], si[[ycol]])))
    mul <- r * 0.8 / max(sqrt(rowSums(bp[, 1:2]^2)))
    bp$x <- bp[[1]] * mul; bp$y <- bp[[2]] * mul
    bp$name <- rownames(bp)
    p <- p +
      ggplot2::geom_segment(data = bp,
                            ggplot2::aes(x = 0, y = 0, xend = x, yend = y),
                            colour = "grey25", linewidth = 0.45,
                            arrow = grid::arrow(length = grid::unit(1.6, "mm"),
                                                type = "closed"),
                            inherit.aes = FALSE) +
      ggrepel::geom_text_repel(data = bp, ggplot2::aes(x, y, label = name),
                               size = 2.1, colour = "grey15", seed = 1,
                               inherit.aes = FALSE, show.legend = FALSE)
  }
  # 
  if (show_species && !is.null(sc$species) && nrow(sc$species) >= 1) {
    sp <- as.data.frame(sc$species)
    sp$name <- rownames(sp)
    sp <- sp[order(rowSums(sp[, 1:2, drop = FALSE]^2), decreasing = TRUE), ]
    sp <- utils::head(sp, 6)
    sp$x <- sp[[1]]; sp$y <- sp[[2]]
    p <- p +
      ggplot2::geom_point(data = sp, ggplot2::aes(x, y),
                          shape = 3, size = 1.1, colour = "grey35",
                          inherit.aes = FALSE, show.legend = FALSE) +
      ggrepel::geom_text_repel(data = sp,
                               ggplot2::aes(x, y, label = name),
                               size = 1.8, colour = "grey25", seed = 1,
                               inherit.aes = FALSE, show.legend = FALSE)
  }
  attr(p, "model") <- ord
  p
}

#' Environmental vector fitting on an existing PCoA
#'
#'  PCoA  `vegan::envfit` 
#'  R2 ""
#'
#' @inheritParams plot_ordination_constrained
#' @param dist Beta distance key from NosZRef output, e.g. `"bray_curtis"`.
#' @return A ggplot object.
#' @export
plot_envfit_ordination <- function(x, branch, env, dist = NULL,
                                   palette = "nosz", permutations = 999) {
  stopifnot(inherits(x, "noszref"))
  env <- prepare_env(x, env)
  if (is.null(env)) nosz_abort("")
  beta <- nosz_need(x, "beta", branch)
  if (is.null(dist)) dist <- grep("bray", names(beta), value = TRUE)[1] %||% names(beta)[1]
  m <- beta[[dist]]
  if (is.null(m)) nosz_abort(" `", dist, "` ")
  smp <- intersect(rownames(m), env$SampleID)
  m <- m[smp, smp, drop = FALSE]
  env <- env[match(smp, env$SampleID), , drop = FALSE]
  rownames(env) <- env$SampleID
  env$SampleID <- NULL
  pc <- stats::cmdscale(stats::as.dist(m), k = 2, eig = TRUE)
  expl <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
  sc <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2],
                   SampleID = smp, stringsAsFactors = FALSE)
  sc <- merge(sc, x$design[, c("SampleID", "Group")], by = "SampleID", all.x = TRUE)
  ef <- try(vegan::envfit(pc$points, env, permutations = permutations), silent = TRUE)
  cols <- group_palette(sc$Group, palette)

  p <- ggplot2::ggplot(sc, ggplot2::aes(PC1, PC2, colour = Group)) +
    ggplot2::geom_point(size = 2, alpha = 0.95) +
    ggplot2::stat_ellipse(ggplot2::aes(fill = Group), geom = "polygon",
                          level = 0.68, alpha = 0.12, colour = NA,
                          show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = sprintf("PCoA1 (%.1f%%)", expl[1]),
                  y = sprintf("PCoA2 (%.1f%%)", expl[2]),
                  title = paste0(branch, " PCoA with environmental fit (", dist, ")"),
                  subtitle = paste0("arrows: R2 / P shown for P < 0.1")) +
    theme_nosz()
  if (!inherits(ef, "try-error")) {
    v <- ef$vectors$arrows[, 1:2, drop = FALSE]
    r2 <- ef$vectors$r
    pv <- ef$vectors$pvals
    df <- data.frame(v, r2 = r2, p = pv, name = rownames(v))
    df <- df[df$p < 0.1, , drop = FALSE]
    if (nrow(df) > 0) {
      rg <- max(abs(c(sc$PC1, sc$PC2)))
      mul <- rg * 0.85 / max(sqrt(rowSums(df[, 1:2]^2)))
      df$AR1 <- df[[1]] * mul; df$AR2 <- df[[2]] * mul
      df$lab <- sprintf("%s (R2=%.2f, %s)", df$name, df$r2, p_star(df$p))
      p <- p +
        ggplot2::geom_segment(data = df,
                              ggplot2::aes(x = 0, y = 0, xend = AR1, yend = AR2),
                              colour = "grey20", linewidth = 0.45,
                              arrow = grid::arrow(length = grid::unit(1.6, "mm"),
                                                  type = "closed"),
                              inherit.aes = FALSE) +
        ggrepel::geom_text_repel(data = df, ggplot2::aes(AR1, AR2, label = lab),
                                 size = 1.9, colour = "grey10", seed = 1,
                                 inherit.aes = FALSE, show.legend = FALSE)
    }
  }
  p
}

#'
# LEfSe 
#'

#' LEfSe-style differential analysis
#'
#'  LEfSe Kruskal-Wallis  ->  Wilcoxon  ->
#'  log10 
#'  LEfSe 
#' "LEfSe-style analysis implemented in nosZviz"
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param kw_p Raw P cut for the Kruskal-Wallis screening step.
#' @param wilcox_p Raw P cut for the pairwise Wilcoxon step.
#' @param lda_cut Absolute effect-size cut (log10 mean-ratio).
#' @return A data frame with `Taxon`, `Group` (winning group), `LDA` (signed),
#'   `mean_in_group`, `mean_rest`, `kw_p`, `wilcox_p`.
#' @export
calc_lefse <- function(x, branch, level = "genus", kw_p = 0.05,
                       wilcox_p = 0.05, lda_cut = 1) {
  stopifnot(inherits(x, "noszref"))
  m <- nosz_tax_matrix(x, branch, level, value = "rpkm")
  m <- rel_abund(m)
  d <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
  m <- m[, d$SampleID, drop = FALSE]
  m <- m[rowMeans(m > 0) >= 0.2, , drop = FALSE]
  gs <- unique(d$Group)
  if (length(gs) < 2) nosz_abort("LEfSe ")
  out <- list()
  for (i in seq_len(nrow(m))) {
    v <- m[i, ]
    kw <- suppressWarnings(stats::kruskal.test(v, factor(d$Group)))
    if (is.na(kw$p.value) || kw$p.value >= kw_p) next
    gm <- tapply(v, d$Group, mean)
    #  vs 
    best <- names(sort(gm, decreasing = TRUE))[1]
    rest <- setdiff(gs, best)
    ok <- all(vapply(rest, function(g) {
      w <- suppressWarnings(stats::wilcox.test(v[d$Group == best],
                                               v[d$Group == g], exact = FALSE))
      !is.na(w$p.value) && w$p.value < wilcox_p
    }, logical(1)))
    if (!ok) next
    lda <- log10((mean(v[d$Group == best]) + 1e-7) /
                   (mean(v[d$Group != best]) + 1e-7))
    if (abs(lda) < lda_cut) next
    out[[length(out) + 1]] <- data.frame(
      Taxon = rownames(m)[i], Group = best, LDA = lda,
      mean_in_group = mean(v[d$Group == best]),
      mean_rest = mean(v[d$Group != best]),
      kw_p = kw$p.value, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) {
    return(data.frame(Taxon = character(), Group = character(),
                      LDA = numeric(), mean_in_group = numeric(),
                      mean_rest = numeric(), kw_p = numeric()))
  }
  do.call(rbind, out)
}

#' LEfSe-style effect size bar plot
#'
#' 
#' LEfSe Nature  nosZ 
#'
#' @param x A `noszref` object.
#' @param res Optional result of [calc_lefse()].
#' @param ... Passed to [calc_lefse()].
#' @return A ggplot object.
#' @export
plot_lefse <- function(x, res = NULL, ...) {
  if (is.null(res)) res <- calc_lefse(x, ...)
  if (nrow(res) == 0) {
    nosz_abort(" LEfSe  kw_p / wilcox_p / lda_cut")
  }
  res <- res[order(res$LDA), , drop = FALSE]
  res$Taxon <- factor(res$Taxon, levels = res$Taxon)
  cols <- group_palette(res$Group, "nosz")
  ggplot2::ggplot(res, ggplot2::aes(Taxon, LDA, fill = Group)) +
    ggplot2::geom_col(width = 0.7, colour = "grey25", linewidth = 0.15) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL, y = "LDA score (log10 mean ratio)",
                  title = "LEfSe-style differential taxa",
                  subtitle = "KW + pairwise Wilcoxon; bar colour = enriched group") +
    theme_nosz(legend_position = "top")
}

#'
# IndVal + 
#'

#' Indicator taxa analysis (IndVal with permutations)
#'
#'  IndVal  A   B+ 
#' ""
#' " biosignature"
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param nperm Number of label permutations.
#' @param value `"rpkm"` or `"relative"`.
#' @return A data frame with `Taxon`, `Group`, `A` (specificity),
#'   `B` (fidelity), `IndVal`, `p`, `padj`.
#' @export
calc_indicator <- function(x, branch, level = "genus", nperm = 999,
                           value = c("rpkm", "relative")) {
  stopifnot(inherits(x, "noszref"))
  value <- match.arg(value)
  m <- nosz_tax_matrix(x, branch, level, value = "rpkm")
  m <- rel_abund(m)
  d <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
  m <- m[, d$SampleID, drop = FALSE]
  m <- m[rowSums(m) > 0, , drop = FALSE]
  gs <- unique(d$Group)
  if (length(gs) < 2) nosz_abort("")

  indval_one <- function(mm, grp) {
    # A:  / B: 
    A <- sapply(gs, function(g) rowMeans(mm[, grp == g, drop = FALSE]))
    A <- A / rowSums(A)
    B <- sapply(gs, function(g) rowMeans(mm[, grp == g, drop = FALSE] > 0))
    A * B
  }
  grp0 <- d$Group
  obs <- indval_one(m, grp0)
  # 
  n <- length(grp0)
  n1 <- sum(grp0 == gs[1])
  if (length(gs) == 2 && choose(n, n1) <= 5000) {
    subsets <- utils::combn(n, n1, simplify = FALSE)
    stats_perm <- lapply(subsets, function(idx) {
      lab <- rep(gs[2], n); lab[idx] <- gs[1]
      indval_one(m, lab)
    })
    cnt <- Reduce(`+`, lapply(stats_perm, function(sp) sp >= obs - 1e-12))
    pv <- cnt / length(subsets)
  } else {
    cnt <- matrix(0, nrow = nrow(obs), ncol = ncol(obs))
    for (i in seq_len(nperm)) {
      prm <- sample(grp0)
      cnt <- cnt + (indval_one(m, prm) >= obs - 1e-12)
    }
    pv <- (cnt + 1) / (nperm + 1)
  }
  # A / B
  Am <- sapply(gs, function(g) rowMeans(m[, grp0 == g, drop = FALSE]))
  Am <- Am / rowSums(Am)
  Bm <- sapply(gs, function(g) rowMeans(m[, grp0 == g, drop = FALSE] > 0))
  out <- data.frame(
    Taxon = rep(rownames(obs), times = ncol(obs)),
    Group = rep(colnames(obs), each = nrow(obs)),
    A = as.numeric(Am),
    B = as.numeric(Bm),
    IndVal = as.numeric(obs),
    p = as.numeric(pv),
    stringsAsFactors = FALSE)
  out$padj <- stats::p.adjust(out$p, "BH")
  out <- out[order(-out$IndVal), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Indicator taxa lollipop plot
#'
#' @param x A `noszref` object.
#' @param res Optional result of [calc_indicator()].
#' @param padj Cut for displayed taxa.
#' @param min_indval Minimum IndVal displayed.
#' @param top_n Maximum number of taxa displayed.
#' @param ... Passed to [calc_indicator()].
#' @return A ggplot object.
#' @export
plot_indicator <- function(x, res = NULL, padj = 0.05, min_indval = 0.25,
                           top_n = 20, ...) {
  if (is.null(res)) res <- calc_indicator(x, ...)
  sig <- res[res$padj < padj & res$IndVal >= min_indval, , drop = FALSE]
  if (nrow(sig) == 0) {
    #  IndVal 
    sig <- utils::head(res[order(-res$IndVal), ], top_n)
    sig <- sig[sig$IndVal > 0, , drop = FALSE]
    if (nrow(sig) == 0) nosz_abort("")
    note <- "no significant indicators at given cut; showing top IndVal taxa"
  } else {
    sig <- utils::head(sig[order(-sig$IndVal), ], top_n)
    note <- sprintf("padj < %g, IndVal >= %g; dot size = fidelity", padj, min_indval)
  }
  res <- sig
  res <- res[order(res$Group, res$IndVal), , drop = FALSE]
  res$Taxon <- factor(res$Taxon, levels = unique(res$Taxon))
  cols <- group_palette(res$Group, "nosz")
  ggplot2::ggplot(res, ggplot2::aes(Taxon, IndVal)) +
    ggplot2::geom_segment(ggplot2::aes(x = Taxon, xend = Taxon, y = 0, yend = IndVal,
                                       colour = Group), linewidth = 0.5,
                          show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(fill = Group, size = B), shape = 21,
                        colour = "grey25", stroke = 0.3) +
    ggplot2::scale_fill_manual(values = cols, name = "Group") +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_size_continuous(range = c(1.2, 3), name = "Fidelity B") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "IndVal (A x B)",
                  title = "Indicator taxa (IndVal)",
                  subtitle = note) +
    theme_nosz(legend_position = "right")
}

#'
#  OTU / 
#'

#' Shared OTU / taxa across groups (UpSet)
#'
#'  OTU/UpSet 
#' " nosZ "
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level `"otu"`  OTU
#' @return An UpSetR plot object.
#' @export
plot_group_upset <- function(x, branch, level = "otu") {
  stopifnot(inherits(x, "noszref"))
  m <- if (identical(level, "otu")) {
    nosz_otu_matrix(x, branch)
  } else {
    nosz_tax_matrix(x, branch, level, value = "count")
  }
  d <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
  gs <- unique(d$Group)
  if (length(gs) < 2) nosz_abort("")
  bin <- sapply(gs, function(g) {
    ss <- intersect(d$SampleID[d$Group == g], colnames(m))
    as.integer(rowSums(m[, ss, drop = FALSE] > 0) > 0)
  })
  #  OTU  sapply 
  if (is.null(dim(bin))) bin <- matrix(bin, nrow = nrow(m), ncol = length(gs))
  colnames(bin) <- gs
  bin <- as.data.frame(bin[apply(bin, 1, sum) > 0, , drop = FALSE])
  if (nrow(bin) == 0) nosz_abort("/")
  UpSetR::upset(bin, sets = gs, sets.bar.color = nosz_pal(length(gs), "nosz"),
                mainbar.y.label = "Shared / unique taxa (intersection size)",
                sets.x.label = "Taxa per group",
                text.scale = c(1.3, 1.1, 1, 0.95, 1.1, 0.9),
                order.by = "freq")
}

#'
# alpha x 
#'

#' Alpha diversity vs environmental correlation heatmap
#'
#'  alpha  Spearman 
#'  alpha 
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param env Environmental table (data frame with `SampleID`, or path).
#' @param method `"spearman"` or `"pearson"`.
#' @return A ggplot object.
#' @export
plot_alpha_env_corr <- function(x, branch, env, method = "spearman") {
  stopifnot(inherits(x, "noszref"))
  env <- prepare_env(x, env)
  if (is.null(env)) nosz_abort("")
  a <- nosz_need(x, "alpha", branch)
  a <- a[a$SampleID %in% env$SampleID, , drop = FALSE]
  env <- env[match(a$SampleID, env$SampleID), , drop = FALSE]
  idx <- setdiff(colnames(a), "SampleID")
  evars <- setdiff(colnames(env), "SampleID")
  rows <- list()
  for (iv in idx) {
    for (ev in evars) {
      ct <- suppressWarnings(stats::cor.test(a[[iv]], env[[ev]], method = method))
      rows[[length(rows) + 1]] <- data.frame(Index = iv, Env = ev,
                                             rho = unname(ct$estimate),
                                             p = ct$p.value,
                                             stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, rows)
  d$Index <- factor(d$Index, levels = rev(idx))
  d$star <- p_star(d$p)
  ggplot2::ggplot(d, ggplot2::aes(Env, Index, fill = rho)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = star), size = 2.4, vjust = 0.75) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                  high = "#B2182B", midpoint = 0,
                                  limits = c(-1, 1), name = method) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0(branch, " alpha diversity x environment"),
                  subtitle = "* p<0.05, ** p<0.01, *** p<0.001") +
    theme_nosz(border = FALSE, x_angle = 30)
}

#'
# -
#'

#' Taxa-environment correlation network
#'
#'  Spearman |rho| 
#' 
#' ""
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level for taxa nodes.
#' @param env Environmental table.
#' @param top_n Top taxa included.
#' @param r_cut Minimum |rho| for an edge.
#' @param p_cut Raw P cut for an edge (Spearman; use cautiously with small n).
#' @return A ggraph object.
#' @export
plot_taxa_env_network <- function(x, branch, env, level = "genus", top_n = 20,
                                  r_cut = 0.6, p_cut = 0.05) {
  stopifnot(inherits(x, "noszref"))
  env <- prepare_env(x, env)
  if (is.null(env)) nosz_abort("")
  m <- nosz_tax_matrix(x, branch, level, value = "rpkm")
  m <- rel_abund(m)
  smp <- intersect(colnames(m), env$SampleID)
  m <- m[, smp, drop = FALSE]
  env <- env[match(smp, env$SampleID), , drop = FALSE]
  if (nrow(m) > top_n) m <- m[order(rowMeans(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  evars <- setdiff(colnames(env), "SampleID")

  edges <- list()
  for (i in seq_len(nrow(m))) {
    for (j in seq_along(evars)) {
      ct <- suppressWarnings(stats::cor.test(m[i, ], env[[evars[j]]],
                                             method = "spearman"))
      if (is.na(ct$estimate)) next
      edges[[length(edges) + 1]] <- data.frame(from = rownames(m)[i],
                                               to = evars[j],
                                               rho = unname(ct$estimate),
                                               p = ct$p.value,
                                               stringsAsFactors = FALSE)
    }
  }
  ed <- do.call(rbind, edges)
  ed$padj <- stats::p.adjust(ed$p, "BH")
  ed <- ed[abs(ed$rho) >= r_cut & ed$p < p_cut, , drop = FALSE]
  if (nrow(ed) == 0) nosz_abort(" r_cut / p_cut - r_cut  p_cut")

  # 
  phy_map <- taxa_phylum_map(x, branch, level)
  nodes <- data.frame(
    name = c(rownames(m), evars),
    type = c(rep("Taxon", nrow(m)), rep("Env", length(evars))),
    stringsAsFactors = FALSE)
  nodes$Phylum <- phy_map[nodes$name]
  nodes$Phylum[is.na(nodes$Phylum)] <- "Unknown"
  nodes$Phylum[nodes$type == "Env"] <- "Environment"
  nodes$Size <- c(rowMeans(m) * 100, rep(max(rowMeans(m)) * 100, length(evars)))

  #  rho/p edge aes 
  g <- igraph::graph_from_data_frame(ed[, c("from", "to", "rho", "p")],
                                     directed = FALSE, vertices = nodes)
  set.seed(3)
  pal <- c(nosz_pal(length(unique(na.omit(phy_map))), "tableau"), "grey35")
  names(pal) <- c(sort(unique(na.omit(phy_map))), "Environment")
  if (!all(nodes$Phylum %in% names(pal))) {
    extra <- setdiff(unique(nodes$Phylum), names(pal))
    pal <- c(pal, stats::setNames(nosz_pal(length(extra), "nosz"), extra))
  }
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(width = abs(rho), colour = rho > 0),
                           alpha = 0.55, show.legend = TRUE) +
    ggraph::scale_edge_width(range = c(0.2, 1.1), guide = "none") +
    ggraph::scale_edge_colour_manual(
      values = c(`TRUE` = "#C0392B", `FALSE` = "#2471A3"),
      labels = c(`TRUE` = "positive", `FALSE` = "negative"), name = "correlation") +
    ggraph::geom_node_point(ggplot2::aes(fill = Phylum, size = Size, shape = type),
                            colour = "grey25", stroke = 0.3) +
    ggplot2::scale_shape_manual(values = c(Taxon = 21, Env = 22), name = NULL) +
    ggplot2::scale_fill_manual(values = pal, name = "Phylum") +
    ggplot2::scale_size_continuous(range = c(1.5, 5), guide = "none") +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 1.7,
                           max.overlaps = 20) +
    ggplot2::labs(title = paste0(branch, " ", level, " x environment network"),
                  subtitle = sprintf("Spearman |rho| >= %g, P < %g", r_cut, p_cut)) +
    theme_nosz(border = FALSE)
}

#  ->   BLAST  lineage 
taxa_phylum_map <- function(x, branch, level = "genus") {
  tl <- try(nosz_need(x, "tax_long", branch), silent = TRUE)
  if (inherits(tl, "try-error") || is.null(tl)) return(stats::setNames(rep(NA_character_, 100), NULL)[0])
  ll <- tax_parse(tl$lineage)
  # OTU/ lineage qseqid  SampleID_OTUIDlineage 
  #  lineage  g_ / f_ 
  ph <- ll$Phylum
  ph[is.na(ph)] <- "Unknown"
  #  OTU-> qseqid 
  #  lineage 
  taxa <- rownames(nosz_tax_matrix(x, branch, level, value = "rpkm"))
  # 
  #  genus  lineage 
  gl <- ll$Genus
  out <- vapply(taxa, function(t) {
    hit <- which(!is.na(gl) & gl == t)
    if (length(hit) == 0) {
      # 
      hit <- which(!is.na(gl) & grepl(t, gl, fixed = TRUE))
    }
    if (length(hit) == 0) return(NA_character_)
    names(sort(table(ph[hit]), decreasing = TRUE))[1]
  }, character(1))
  stats::setNames(out, taxa)
}

# UPGMA clustering section
#' Taxa abundance-profile dendrogram
#'
#' Spearman  UPGMA 
#' ""
#'
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param top_n Top taxa shown.
#' @param circular Logical, circular layout.
#' @param label_top Label the top N most abundant leaves.
#' @return A ggraph object.
#' @export
plot_tax_dendrogram <- function(x, branch, level = "genus", top_n = 40,
                                circular = TRUE, label_top = 8) {
  stopifnot(inherits(x, "noszref"))
  m <- nosz_tax_matrix(x, branch, level, value = "rpkm")
  m <- rel_abund(m)
  if (nrow(m) > top_n) m <- m[order(rowMeans(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  if (nrow(m) < 4) nosz_abort("<4")
  cc <- suppressWarnings(stats::cor(t(m), method = "spearman"))
  cc[is.na(cc)] <- 0
  D <- stats::as.dist(1 - cc)
  hc <- stats::hclust(D, method = "average")
  lay <- ggraph::create_layout(hc, layout = "dendrogram", circular = circular)
  lay$Mean <- rowMeans(m)[lay$label]
  phy_map <- taxa_phylum_map(x, branch, level)
  lay$Phylum <- phy_map[lay$label]
  lay$Phylum[is.na(lay$Phylum)] <- "Unknown"
  if (circular) lay$Phylum[!lay$leaf] <- NA
  labs <- lay[lay$leaf, , drop = FALSE]
  labs <- labs[order(-labs$Mean), , drop = FALSE]
  labs <- utils::head(labs, label_top)

  g <- ggraph::ggraph(lay) +
    ggraph::geom_edge_elbow(colour = "grey55", linewidth = 0.3) +
    ggraph::geom_node_point(ggplot2::aes(filter = leaf, size = Mean, fill = Phylum),
                            shape = 21, colour = "grey25", stroke = 0.3,
                            show.legend = TRUE) +
    ggplot2::scale_size_continuous(range = c(1, 4.5), name = "Mean rel. ab.") +
    ggplot2::scale_fill_manual(values = c(nosz_pal(max(length(unique(na.omit(lay$Phylum))) - 1, 1), "tableau"),
                                          "grey85"), na.value = "transparent", name = "Phylum") +
    ggplot2::theme_void(base_size = 8)
  if (nrow(labs) > 0) {
    g <- g + ggraph::geom_node_text(ggplot2::aes(x = x * 1.08, y = y * 1.08,
                                                 filter = label %in% labs$label,
                                                 label = label),
                                    size = 1.7, max.overlaps = 30)
  }
  g + ggplot2::labs(title = paste0(branch, " ", level,
                                   " abundance-profile dendrogram (UPGMA)"),
                    subtitle = "Similarity of abundance profiles across samples - not a phylogeny")
}
