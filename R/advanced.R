# ---------------------------------------------------------------------------
# 模块九：环境因子关联与进阶分析（NosZRef 内置可视化未覆盖）
#  - 约束排序 RDA / CCA / db-RDA（环境因子驱动支系类群分异）
#  - envfit 环境向量叠加（无约束 PCoA 上）
#  - LEfSe 风格差异分析（Kruskal-Wallis + Wilcoxon + 效应量）
#  - 指示种分析（IndVal + 置换检验）
#  - 组间共享 OTU UpSet
#  - alpha 指数 × 环境因子相关热图
#  - 类群-环境因子相关网络
#  - 分类学树状图（无树文件时的替代：丰度谱系 UPGMA 环形树）
# ---------------------------------------------------------------------------

# 环境表清洗：匹配样本、仅保留数值列、剔除零方差列
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
#' nosZ 特化的约束排序（RDA / CCA / db-RDA）：以支系内分类单元丰度为
#' 响应矩阵、环境因子为解释变量，封装 vegan 的 rda()/cca()/capscale()，
#' 自动做 Hellinger 转换、剔除零方差变量，并在副标题标注约束轴解释率与
#' 置换检验 P 值。排序图回答：**什么环境条件驱动了该支系类群的分异？**
#'
#' @param x A `noszref` object.
#' @param branch Branch name (one branch at a time).
#' @param level Taxonomic level used as response matrix.
#' @param env Environmental data: data frame with a `SampleID` column and
#'   numeric variables, or path to such a table.
#' @param method `"RDA"`（Hellinger 转换后的冗余分析）, `"CCA"` 或
#'   `"dbRDA"`（基于距离的冗余分析）.
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
  if (is.null(env)) nosz_abort("环境表不可用：需要含 SampleID 列与至少一个非零方差数值列。")
  m <- nosz_tax_matrix(x, branch, level, value = "rpkm")
  smp <- intersect(colnames(m), env$SampleID)
  if (length(smp) < 4) nosz_abort("样本与环境表的重叠少于 4 个，无法做约束排序。")
  m <- m[, smp, drop = FALSE]
  env <- env[match(smp, env$SampleID), , drop = FALSE]
  rownames(env) <- env$SampleID
  env$SampleID <- NULL
  # top_n 分类单元进入模型；转置为 样本 x 类群（vegan 要求行=样点）
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
  # 残差秩为 0 时可能只有 1 个轴：补第二轴为 0 以便绘图
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

  # 环境因子箭头（双序图分数，缩放到数据范围 0.8；vegan 返回键名为 biplot）
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
  # 高载荷类群标签
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
#' 在无约束 PCoA 上叠加 `vegan::envfit` 拟合的环境因子向量，
#' 副标题标注每个变量的 R2 与显著性，回答"群落分异沿哪些环境梯度展开"。
#'
#' @inheritParams plot_ordination_constrained
#' @param dist Beta distance key from NosZRef output, e.g. `"bray_curtis"`.
#' @return A ggplot object.
#' @export
plot_envfit_ordination <- function(x, branch, env, dist = NULL,
                                   palette = "nosz", permutations = 999) {
  stopifnot(inherits(x, "noszref"))
  env <- prepare_env(x, env)
  if (is.null(env)) nosz_abort("环境表不可用。")
  beta <- nosz_need(x, "beta", branch)
  if (is.null(dist)) dist <- grep("bray", names(beta), value = TRUE)[1] %||% names(beta)[1]
  m <- beta[[dist]]
  if (is.null(m)) nosz_abort("距离矩阵 `", dist, "` 不存在。")
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

# ---------------------------------------------------------------------------
# LEfSe 风格差异分析
# ---------------------------------------------------------------------------

#' LEfSe-style differential analysis
#'
#' 简化版 LEfSe 流程（Kruskal-Wallis 组间筛选 -> 组内两两 Wilcoxon 精简 ->
#' 效应量以两组相对丰度均值之比的 log10 计）。效应量方向指向均值最高的
#' 分组。注意：这是 LEfSe 的近似实现，投稿时请在方法中写明
#' "LEfSe-style analysis implemented in nosZviz"。
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
  if (length(gs) < 2) nosz_abort("LEfSe 至少需要两个分组。")
  out <- list()
  for (i in seq_len(nrow(m))) {
    v <- m[i, ]
    kw <- suppressWarnings(stats::kruskal.test(v, factor(d$Group)))
    if (is.na(kw$p.value) || kw$p.value >= kw_p) next
    gm <- tapply(v, d$Group, mean)
    # 均值最高组 vs 其余合并的两两检验
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
#' 横向条形图展示差异类群的效应量，条形按优势分组着色——
#' LEfSe 的标准呈现方式（Nature 系期刊 nosZ 论文常用）。
#'
#' @param x A `noszref` object.
#' @param res Optional result of [calc_lefse()].
#' @param ... Passed to [calc_lefse()].
#' @return A ggplot object.
#' @export
plot_lefse <- function(x, res = NULL, ...) {
  if (is.null(res)) res <- calc_lefse(x, ...)
  if (nrow(res) == 0) {
    nosz_abort("没有通过 LEfSe 阈值的类群；可放宽 kw_p / wilcox_p / lda_cut。")
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

# ---------------------------------------------------------------------------
# 指示种分析（IndVal + 置换检验）
# ---------------------------------------------------------------------------

#' Indicator taxa analysis (IndVal with permutations)
#'
#' 经典 IndVal 指数（特异性 A × 忠实度 B）+ 组标签置换检验，
#' 找出"专属"某分组的高保真类群。比单纯差异检验更能回答
#' "哪个类群可以作为该生境的 biosignature"。
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
  if (length(gs) < 2) nosz_abort("指示种分析至少需要两个分组。")

  indval_one <- function(mm, grp) {
    # A: 组均值 / 各组均值之和；B: 组内出现频率
    A <- sapply(gs, function(g) rowMeans(mm[, grp == g, drop = FALSE]))
    A <- A / rowSums(A)
    B <- sapply(gs, function(g) rowMeans(mm[, grp == g, drop = FALSE] > 0))
    A * B
  }
  grp0 <- d$Group
  obs <- indval_one(m, grp0)
  # 置换检验：两组且排列总数不大时完全枚举（含恒等排列），否则随机置换
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
  # A：特异性（组均值 / 各组均值之和）；B：忠实度（组内出现频率）
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
    # 小样本数据可能无显著指示种：降级展示 IndVal 最高的类群并注明
    sig <- utils::head(res[order(-res$IndVal), ], top_n)
    sig <- sig[sig$IndVal > 0, , drop = FALSE]
    if (nrow(sig) == 0) nosz_abort("指示种表为空。")
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

# ---------------------------------------------------------------------------
# 组间共享 OTU / 类群
# ---------------------------------------------------------------------------

#' Shared OTU / taxa across groups (UpSet)
#'
#' 组间 OTU（或任意层级类群）存在/缺失的集合交叉，UpSet 呈现——
#' 回答"各分组的核心 nosZ 类群与特有类群各有多少"。
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level `"otu"` 或分类层级（用类群存在性代替 OTU）。
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
  if (length(gs) < 2) nosz_abort("至少需要两个分组。")
  bin <- sapply(gs, function(g) {
    ss <- intersect(d$SampleID[d$Group == g], colnames(m))
    as.integer(rowSums(m[, ss, drop = FALSE] > 0) > 0)
  })
  # 单 OTU 时 sapply 退化为向量，兜底成矩阵
  if (is.null(dim(bin))) bin <- matrix(bin, nrow = nrow(m), ncol = length(gs))
  colnames(bin) <- gs
  bin <- as.data.frame(bin[apply(bin, 1, sum) > 0, , drop = FALSE])
  if (nrow(bin) == 0) nosz_abort("没有可比较的存在/缺失模式。")
  UpSetR::upset(bin, sets = gs, sets.bar.color = nosz_pal(length(gs), "nosz"),
                mainbar.y.label = "Shared / unique taxa (intersection size)",
                sets.x.label = "Taxa per group",
                text.scale = c(1.3, 1.1, 1, 0.95, 1.1, 0.9),
                order.by = "freq")
}

# ---------------------------------------------------------------------------
# alpha x 环境因子相关热图
# ---------------------------------------------------------------------------

#' Alpha diversity vs environmental correlation heatmap
#'
#' 六种 alpha 指数与环境因子的 Spearman 相关矩阵（带显著性星号），
#' 是 alpha 多样性与理化性质关联的标准呈现。
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
  if (is.null(env)) nosz_abort("环境表不可用。")
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

# ---------------------------------------------------------------------------
# 类群-环境因子相关网络
# ---------------------------------------------------------------------------

#' Taxa-environment correlation network
#'
#' 类群（行）与环境因子（列）的 Spearman 相关网络：|rho| 与显著性双重
#' 筛边，环境因子为方形节点、类群为圆形节点（按门着色、按丰度定大小），
#' 正负相关用边色区分——快速定位"哪些类群与哪些环境梯度耦合"。
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
  if (is.null(env)) nosz_abort("环境表不可用。")
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
  if (nrow(ed) == 0) nosz_abort("没有通过 r_cut / p_cut 的类群-环境相关；可降低 r_cut 或放宽 p_cut。")

  # 门水平着色映射
  phy_map <- taxa_phylum_map(x, branch, level)
  nodes <- data.frame(
    name = c(rownames(m), evars),
    type = c(rep("Taxon", nrow(m)), rep("Env", length(evars))),
    stringsAsFactors = FALSE)
  nodes$Phylum <- phy_map[nodes$name]
  nodes$Phylum[is.na(nodes$Phylum)] <- "Unknown"
  nodes$Phylum[nodes$type == "Env"] <- "Environment"
  nodes$Size <- c(rowMeans(m) * 100, rep(max(rowMeans(m)) * 100, length(evars)))

  # 边属性 rho/p 必须随图对象一起传入，edge aes 才能取到
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

# 类群 -> 门 映射（从 BLAST 长表 lineage 聚合）
taxa_phylum_map <- function(x, branch, level = "genus") {
  tl <- try(nosz_need(x, "tax_long", branch), silent = TRUE)
  if (inherits(tl, "try-error") || is.null(tl)) return(stats::setNames(rep(NA_character_, 100), NULL)[0])
  ll <- tax_parse(tl$lineage)
  # OTU/类群名与 lineage 的对应：qseqid 形如 SampleID_OTUID，lineage 末位为种
  # 用 lineage 中 g_ / f_ 等名字做名字匹配成本高；这里按频率取最常见门
  ph <- ll$Phylum
  ph[is.na(ph)] <- "Unknown"
  # 长表中无法直接得到 OTU->门 的键（qseqid 与分类矩阵行名不一致时），
  # 退化为按门频率分配：对每个类群名在 lineage 各字段中模糊匹配
  taxa <- rownames(nosz_tax_matrix(x, branch, level, value = "rpkm"))
  # 尝试直接从各层级矩阵建立映射：类群在各门矩阵中的对应关系未知，
  # 采用启发式：该类群名若出现在某门的 genus 名单中（长表 lineage 解析）
  gl <- ll$Genus
  out <- vapply(taxa, function(t) {
    hit <- which(!is.na(gl) & gl == t)
    if (length(hit) == 0) {
      # 模糊：包含关系
      hit <- which(!is.na(gl) & grepl(t, gl, fixed = TRUE))
    }
    if (length(hit) == 0) return(NA_character_)
    names(sort(table(ph[hit]), decreasing = TRUE))[1]
  }, character(1))
  stats::setNames(out, taxa)
}

# ---------------------------------------------------------------------------
# 分类学树状图（丰度谱系 UPGMA，无需树文件的替代方案）
# ---------------------------------------------------------------------------

#' Taxa abundance-profile dendrogram
#'
#' 无系统发育树文件时的替代方案：以类群跨样本丰度谱的相似性
#' （Spearman 相关距离）做 UPGMA 聚类，环形布局；叶节点大小映射平均
#' 相对丰度、颜色映射所属门，可直接读出"同涨同落"的类群模块。
#' 注意：这是**丰度谱系树**而非系统发育树，图注务必写明。
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
  if (nrow(m) < 4) nosz_abort("类群数量过少（<4），无法建树。")
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
