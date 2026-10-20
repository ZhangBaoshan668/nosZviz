# ---------------------------------------------------------------------------
# 模块十：扩展图型（整合"R语言数据分析指南"等公众号/期刊常见可视化案例）
#  - plot_mantel_env()      Mantel test 组合图（环境 Spearman 热图 + Mantel 连线）
#  - plot_chord()           circlize 弦图（分组 <-> 支系丰度流向）
#  - plot_circular_bar()    环形柱状图（分组 x 支系丰度）
#  - plot_circular_heatmap()环形热图（样本 x 类群 z-score）
#  - plot_alpha_ridge()     山脊图（alpha 指数组间分布）
#  - plot_env_corrplot()    环境因子相关 corrplot（经典四色圆圈风格）
#  - plot_taxa_bump()       排名变化 bump 图（类群丰度排名随分组变化）
# ---------------------------------------------------------------------------

#' Mantel test combo plot (environment heatmap + Mantel links)
#'
#' 公众号/高分期刊高频出现的组合图：左侧为环境因子间的 Spearman 相关
#' 下三角热图（带显著性星号），右侧以曲线连接各环境因子与各 nosZ
#' 分支的类群距离矩阵，连线颜色编码 Mantel r，线宽编码 -log10(P)。
#' 回答：**哪些环境因子既彼此相关，又整体上与某支系群落结构耦合？**
#'
#' @param x A `noszref` object.
#' @param env Environmental table (data frame with `SampleID`, or path).
#' @param level Taxonomic level for the community matrices.
#' @param method Distance for community matrices, passed to
#'   [vegan::vegdist()] (on Hellinger-transformed relative abundances).
#' @param permutations Mantel permutations.
#' @param palette Colour palette for Mantel r gradient.
#' @return A patchwork object.
#' @export
plot_mantel_env <- function(x, env, level = "genus",
                            method = c("bray", "jaccard", "euclidean"),
                            permutations = 999,
                            palette = c("#2166AC", "white", "#B2182B")) {
  stopifnot(inherits(x, "noszref"))
  method <- match.arg(method)
  env <- prepare_env(x, env)
  if (is.null(env)) nosz_abort("环境表不可用。")
  evars <- setdiff(colnames(env), "SampleID")
  if (length(evars) < 2) nosz_abort("环境表至少需要两个数值变量。")

  # ---- 左：环境间 Spearman 下三角热图 ----
  em <- env[, evars, drop = FALSE]
  cors <- list()
  for (i in seq_along(evars)) for (j in seq_along(evars)) {
    if (i >= j) next
    ct <- suppressWarnings(stats::cor.test(em[[i]], em[[j]], method = "spearman"))
    cors[[length(cors) + 1]] <- data.frame(Var1 = evars[i], Var2 = evars[j],
                                           rho = unname(ct$estimate),
                                           p = ct$p.value,
                                           stringsAsFactors = FALSE)
  }
  cd <- do.call(rbind, cors)
  cd$star <- p_star(cd$p)
  lv <- unique(c(cd$Var1, cd$Var2))
  cd$Var1 <- factor(cd$Var1, levels = rev(lv))
  cd$Var2 <- factor(cd$Var2, levels = lv)

  p_heat <- ggplot2::ggplot(cd, ggplot2::aes(Var2, Var1, fill = rho)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = star), size = 2.6, vjust = 0.8) +
    ggplot2::scale_fill_gradient2(low = palette[1], mid = palette[2],
                                  high = palette[3], midpoint = 0,
                                  limits = c(-1, 1), name = "Spearman rho") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_nosz(border = FALSE, x_angle = 45) +
    ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 14, 4, "pt"))

  # ---- 右：Mantel 连线 ----
  mant <- list()
  for (b in x$branches) {
    m <- try(nosz_tax_matrix(x, b, level, value = "rpkm"), silent = TRUE)
    if (inherits(m, "try-error") || is.null(m)) next
    m <- rel_abund(m)
    smp <- intersect(colnames(m), env$SampleID)
    if (length(smp) < 4) next
    m <- t(m[, smp, drop = FALSE])
    dm <- vegan::vegdist(vegan::decostand(m, "hellinger"), method = method)
    for (ev in evars) {
      de <- stats::dist(env[[ev]][match(smp, env$SampleID)])
      mt <- try(vegan::mantel(dm, de, permutations = permutations), silent = TRUE)
      if (inherits(mt, "try-error")) next
      mant[[length(mant) + 1]] <- data.frame(Branch = b, Env = ev,
                                             r = unname(mt$statistic),
                                             p = mt$signif,
                                             stringsAsFactors = FALSE)
    }
  }
  md <- do.call(rbind, mant)
  if (is.null(md) || nrow(md) == 0) {
    nosz_abort("Mantel 检验无法完成（分支距离矩阵或环境数据不足）。")
  }
  # 布局：环境变量在左侧一列（y 位置与热图行对齐），分支在右侧一列
  n_ev <- length(evars)
  n_br <- length(unique(md$Branch))
  ev_y <- stats::setNames(seq(n_ev, 1), evars)          # 与热图行序一致
  br_y <- stats::setNames(seq(n_br, 1), unique(md$Branch))
  md$y  <- ev_y[md$Env]
  md$yend <- br_y[md$Branch]
  md$x <- 0; md$xend <- 1
  md$nlp <- -log10(pmax(md$p, 1 / (permutations + 1)))

  lab_ev <- data.frame(Env = evars, y = ev_y[evars], x = -0.02)
  lab_br <- data.frame(Branch = unique(md$Branch), y = br_y[unique(md$Branch)],
                       x = 1.02)
  p_link <- ggplot2::ggplot(md, ggplot2::aes(x, y, xend = xend, yend = yend)) +
    ggplot2::geom_curve(ggplot2::aes(colour = r, linewidth = nlp),
                        curvature = -0.35, show.legend = TRUE) +
    ggplot2::geom_text(data = lab_ev, ggplot2::aes(x, y, label = Env),
                       hjust = 1, size = 2.4, inherit.aes = FALSE) +
    ggplot2::geom_text(data = lab_br, ggplot2::aes(x, y, label = Branch),
                       hjust = 0, size = 2.4, fontface = "bold",
                       colour = clade_cols(unique(md$Branch))[unique(md$Branch)],
                       inherit.aes = FALSE) +
    ggplot2::scale_colour_gradient2(low = palette[1], mid = palette[2],
                                    high = palette[3], midpoint = 0,
                                    limits = c(-1, 1), name = "Mantel r") +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 1.3),
                                        name = expression(-log[10] * " P")) +
    ggplot2::coord_cartesian(xlim = c(0, 1), clip = "off") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Mantel: clade community x env") +
    theme_nosz(border = FALSE, legend_position = "right") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(4, 46, 4, 34, "pt"))

  patchwork::wrap_plots(p_heat, p_link, ncol = 2, widths = c(1, 1.15))
}

#' Chord diagram of group-clade abundance flows
#'
#' circlize 弦图：分组（扇区）与 nosZ 分支之间的丰度"流向"，每条弦的
#' 宽度正比于该分组下该支系的平均 RPKM 占比——一眼读出"哪群样本把
#' 丰度贡献给了哪个支系"。
#'
#' @param x A `noszref` object.
#' @param value `"rpkm"` 或 `"relative"`（分组内占比）。
#' @return An object of class `nosz_chord`（print 时即时绘制弦图）.
#' @export
plot_chord <- function(x, value = c("rpkm", "relative")) {
  stopifnot(inherits(x, "noszref"))
  value <- match.arg(value)
  mat <- sapply(x$branches, function(b) {
    r <- x$rpkm[[b]]
    if (is.null(r)) return(rep(NA_real_, nrow(x$design)))
    vapply(x$design$SampleID, function(s) {
      i <- match(s, r$SampleID); if (is.na(i)) 0 else r$RPKM[i]
    }, numeric(1))
  })
  colnames(mat) <- x$branches
  # 分组均值
  gs <- unique(x$design$Group)
  mat <- t(sapply(gs, function(g) colMeans(mat[x$design$Group == g, , drop = FALSE])))
  rownames(mat) <- gs
  if (value == "relative") {
    tot <- sum(mat); if (tot > 0) mat <- mat / tot * 100
  }
  if (sum(mat) == 0) nosz_abort("分支 RPKM 表缺失。")

  grid_cols <- c(stats::setNames(nosz_pal(length(gs), "nosz"), gs),
                 clade_cols(x$branches))
  # circlize/corrplot 为即时绘制图形：返回带 print 方法的对象，
  # 使其既能在交互窗口直接 print，也能被 nosz_save() 正确输出
  structure(list(mat = mat, cols = grid_cols), class = "nosz_chord")
}

#' @export
print.nosz_chord <- function(x, ...) {
  # 出错也要清理 circos 状态，避免污染后续绘图
  on.exit(circlize::circos.clear(), add = TRUE)
  circlize::chordDiagram(x$mat, grid.col = x$cols,
                         annotationTrack = "grid", preAllocateTracks = 1,
                         transparency = 0.35)
  circlize::circos.trackPlotRegion(track.index = 1,
                                   panel.fun = function(x, y) {
    # panel.fun 的形参是数据坐标；扇区名从 cell 元数据获取
    xlim <- circlize::get.cell.meta.data("xlim")
    ylim <- circlize::get.cell.meta.data("ylim")
    sector.name <- circlize::get.cell.meta.data("sector.index")
    circlize::circos.text(mean(xlim), ylim[2] + 0.12, sector.name,
                          facing = "clockwise", niceFacing = TRUE,
                          adj = c(0, 0.5), cex = 0.55)
  }, bg.border = NA)
  invisible(x)
}

#' Circular bar plot of group-level clade abundance
#'
#' 环形柱状图：外圈按分组分面、支系着色的组均值 RPKM/占比柱，
#' 公众号常见的"高颜值组成图"变体。
#'
#' @param x A `noszref` object.
#' @param value `"rpkm"` 或 `"relative"`。
#' @param label_size Label font size.
#' @return A ggplot object.
#' @export
plot_circular_bar <- function(x, value = c("rpkm", "relative"),
                              label_size = 2.4) {
  stopifnot(inherits(x, "noszref"))
  value <- match.arg(value)
  d <- calc_clade_stats(x)
  gs <- unique(d$Group)
  long <- do.call(rbind, lapply(gs, function(g) {
    dd <- d[d$Group == g, , drop = FALSE]
    tot <- sum(dd$Total_RPKM, na.rm = TRUE)
    do.call(rbind, lapply(x$branches, function(b) {
      v <- sum(dd[[paste0(b, "_RPKM")]], na.rm = TRUE)
      data.frame(Group = g, Branch = b,
                 Value = if (value == "relative") (if (tot > 0) v / tot * 100 else 0) else v,
                 stringsAsFactors = FALSE)
    }))
  }))
  if (sum(long$Value) == 0) nosz_abort("分支 RPKM 表缺失。")
  cols <- clade_cols(x$branches)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = Group, y = Value, fill = Branch)) +
    ggplot2::geom_col(width = 0.78, colour = "grey20", linewidth = 0.15,
                      position = ggplot2::position_dodge(width = 0.85)) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::coord_polar() +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0("nosZ clade abundance by group (", value, ")"),
                  subtitle = "Bars dodge within each group sector") +
    theme_nosz(border = FALSE) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey92",
                                                              linewidth = 0.2))
  p
}

#' Circular heatmap of sample-by-taxon profiles
#'
#' 环形热图：样本为列、类群为行（z-score），按支系分扇区排布，
#' 中心留白成环。适合作为"多分支一图流"的展示页首图。
#'
#' @param x A `noszref` object.
#' @param level Taxonomic level.
#' @param top_n Top taxa per branch.
#' @param palette Diverging palette (low/mid/high) for z-scores.
#' @return A ggplot object.
#' @export
plot_circular_heatmap <- function(x, level = "genus", top_n = 12,
                                  palette = c("#2166AC", "white", "#B2182B")) {
  stopifnot(inherits(x, "noszref"))
  blocks <- list()
  for (b in x$branches) {
    m <- try(nosz_tax_matrix(x, b, level, value = "rpkm"), silent = TRUE)
    if (inherits(m, "try-error") || is.null(m)) next
    m <- rel_abund(m)
    if (nrow(m) > top_n) m <- m[order(rowMeans(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
    z <- zscore_rows(m)
    long <- data.frame(Branch = b, Taxon = rep(rownames(z), times = ncol(z)),
                       SampleID = rep(colnames(z), each = nrow(z)),
                       z = as.numeric(z), stringsAsFactors = FALSE)
    blocks[[length(blocks) + 1]] <- long
  }
  d <- do.call(rbind, blocks)
  if (is.null(d) || nrow(d) == 0) nosz_abort("没有可用的 ", level, " 层级矩阵。")
  d$Branch <- factor(d$Branch, levels = x$branches)
  d$Taxon <- factor(d$Taxon, levels = rev(unique(d$Taxon)))
  smp_ord <- order_samples(unique(d$SampleID), x$design)
  d$SampleID <- factor(d$SampleID, levels = smp_ord)
  lim <- max(2, stats::quantile(abs(d$z), 0.98, na.rm = TRUE))

  ggplot2::ggplot(d, ggplot2::aes(SampleID, Taxon, fill = z)) +
    ggplot2::geom_tile(colour = NA) +
    ggplot2::scale_fill_gradientn(
      colours = colorRampPalette(palette)(101), limits = c(-lim, lim),
      name = "z-score") +
    ggplot2::coord_polar() +
    ggplot2::facet_wrap(~Branch, ncol = 1) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0("Sample x ", level,
                                 " circular heatmap (z-score, by clade)")) +
    theme_nosz(border = FALSE) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5),
                   axis.text.x = ggplot2::element_text(size = 6),
                   strip.text = ggplot2::element_text(face = "bold"))
}

#' Ridgeline plot of alpha diversity
#'
#' 山脊图（ggridges）：alpha 指数按分组的密度脊线堆叠，
#' 直观展示组间分布的偏移与拖尾——公众号常见的分布比较图型。
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param index Alpha index, e.g. `"shannon"`.
#' @param palette Palette for groups.
#' @return A ggplot object.
#' @export
plot_alpha_ridge <- function(x, branch, index = "shannon", palette = "nosz") {
  stopifnot(inherits(x, "noszref"))
  a <- nosz_need(x, "alpha", branch)
  if (!index %in% colnames(a)) {
    nosz_abort("指数 `", index, "` 不在 alpha 表中。可用：",
               paste(setdiff(colnames(a), "SampleID"), collapse = ", "))
  }
  d <- merge(a[, c("SampleID", index)], x$design[, c("SampleID", "Group")],
             by = "SampleID", all.x = TRUE)
  d$Group <- factor(d$Group, levels = rev(sort(unique(d$Group))))
  cols <- group_palette(d$Group, palette)
  ggplot2::ggplot(d, ggplot2::aes(.data[[index]], Group, fill = Group)) +
    ggridges::geom_density_ridges(rel_min_height = 0.01, alpha = 0.75,
                                  colour = "grey25",
                                  scale = 1.05) +
    ggplot2::geom_point(ggplot2::aes(colour = Group), shape = 124,
                        size = 1.6, stroke = 1.4, alpha = 0.6,
                        position = ggplot2::position_jitter(width = 0,
                                                            height = 0.25,
                                                            seed = 1),
                        show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(x = index, y = NULL,
                  title = paste0(branch, " ", index, " ridgeline by group")) +
    theme_nosz()
}

#' Environment correlation heatmap (ggplot)
#'
#' 环境因子间 Spearman 相关矩阵热图（上三角 + 数值 + 显著性星号），
#' ggplot 实现保证与包内其它图型风格统一、跨设备稳定。
#'
#' @param x A `noszref` object.
#' @param env Environmental table.
#' @param method Correlation method.
#' @return A ggplot object.
#' @export
plot_env_corrplot <- function(x, env, method = "spearman") {
  stopifnot(inherits(x, "noszref"))
  env <- prepare_env(x, env)
  if (is.null(env)) nosz_abort("环境表不可用。")
  evars <- setdiff(colnames(env), "SampleID")
  if (length(evars) < 2) nosz_abort("环境表至少需要两个数值变量。")
  n <- length(evars)
  cells <- expand.grid(Var1 = evars, Var2 = evars, stringsAsFactors = FALSE)
  # 上三角（含对角线）：Var2 的位置序号 >= Var1
  cells <- cells[match(cells$Var2, evars) >= match(cells$Var1, evars), , drop = FALSE]
  cells$rho <- NA_real_; cells$p <- NA_real_
  for (i in seq_len(nrow(cells))) {
    a <- cells$Var1[i]; b <- cells$Var2[i]
    if (identical(a, b)) { cells$rho[i] <- 1; cells$p[i] <- NA; next }
    ct <- suppressWarnings(stats::cor.test(env[[a]], env[[b]], method = method))
    cells$rho[i] <- unname(ct$estimate)
    cells$p[i] <- ct$p.value
  }
  cells$lab <- ifelse(is.na(cells$rho), "",
                      sprintf("%.2f", cells$rho))
  cells$star <- p_star(cells$p)
  cells$lab <- paste0(cells$lab, ifelse(nzchar(cells$star),
                                        paste0(" ", cells$star), ""))
  cells$Var1 <- factor(cells$Var1, levels = rev(evars))
  cells$Var2 <- factor(cells$Var2, levels = evars)
  ggplot2::ggplot(cells, ggplot2::aes(Var2, Var1, fill = rho)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 2.2,
                       colour = "grey15") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                  high = "#B2182B", midpoint = 0,
                                  limits = c(-1, 1), name = method) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Environmental variable correlation",
                  subtitle = "* p<0.05, ** p<0.01, *** p<0.001") +
    theme_nosz(border = FALSE, x_angle = 30)
}

#' Bump chart of taxa rank changes across groups
#'
#' 排名变化图（bump chart）：Top 类群在各分组中的平均丰度排名随分组
#' 排布的折线——展示"谁在哪个分组里登顶、谁跌落"。
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param level Taxonomic level.
#' @param top_n Top taxa tracked.
#' @param palette Palette for taxa lines.
#' @return A ggplot object.
#' @export
plot_taxa_bump <- function(x, branch, level = "genus", top_n = 10,
                           palette = "tableau") {
  stopifnot(inherits(x, "noszref"))
  m <- nosz_tax_matrix(x, branch, level, value = "rpkm")
  m <- rel_abund(m)
  d <- x$design[x$design$SampleID %in% colnames(m), , drop = FALSE]
  gs <- unique(d$Group)
  if (length(gs) < 2) nosz_abort("bump 图至少需要两个分组。")
  gm <- sapply(gs, function(g) {
    ss <- intersect(d$SampleID[d$Group == g], colnames(m))
    rowMeans(m[, ss, drop = FALSE])
  })
  colnames(gm) <- gs
  top <- names(sort(rowMeans(gm), decreasing = TRUE))[seq_len(min(top_n, nrow(gm)))]
  gm <- gm[top, , drop = FALSE]
  rk <- apply(-gm, 2, rank)
  long <- data.frame(Taxon = rep(rownames(rk), times = ncol(rk)),
                     Group = rep(colnames(rk), each = nrow(rk)),
                     Rank = as.numeric(rk), stringsAsFactors = FALSE)
  cols <- stats::setNames(nosz_pal(length(top), palette), top)
  long$Group <- factor(long$Group, levels = gs)
  long$Taxon <- factor(long$Taxon, levels = top)
  ggplot2::ggplot(long, ggplot2::aes(Group, Rank, group = Taxon,
                                     colour = Taxon)) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_line(linewidth = 0.55, alpha = 0.85) +
    ggrepel::geom_text_repel(data = subset(long, Group == gs[1]),
                             ggplot2::aes(label = Taxon), hjust = 1.1,
                             size = 1.7, direction = "y", seed = 1,
                             show.legend = FALSE) +
    ggrepel::geom_text_repel(data = subset(long, Group == gs[length(gs)]),
                             ggplot2::aes(label = Taxon), hjust = -0.1,
                             size = 1.7, direction = "y", seed = 1,
                             show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_y_reverse(breaks = seq_len(nrow(rk))) +
    ggplot2::labs(x = NULL, y = "Rank (by mean relative abundance)",
                  title = paste0(branch, " ", level,
                                 " rank changes across groups")) +
    theme_nosz() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold"))
}
