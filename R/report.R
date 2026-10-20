# ---------------------------------------------------------------------------
# 模块七：综合面板与一键出图
# ---------------------------------------------------------------------------

# 取 OTU 数最多的分支作为默认展示分支
.primary_branch <- function(x, prefer = "nosZI") {
  if (prefer %in% x$branches) return(prefer)
  x$branches[1]
}

#' Overview dashboard
#'
#' A four-panel summary of a NosZRef run: branch abundance, branch partition,
#' alpha diversity and beta diversity ordination.
#'
#' @param x A `noszref` object.
#' @param branch Branch used for the alpha and beta panels.
#' @param palettes Named list of palettes for `abundance`, `composition`,
#'   `alpha` and `beta` panels.
#' @return A patchwork object.
#' @export
plot_overview_dashboard <- function(x, branch = NULL,
                                    palettes = list(abundance = "nosz",
                                                    composition = "branch",
                                                    alpha = "nosz",
                                                    beta = "nosz")) {
  stopifnot(inherits(x, "noszref"))
  branch <- branch %||% .primary_branch(x)
  p1 <- plot_branch_rpkm(x, style = "bar", facet = TRUE,
                         palette = palettes$abundance)
  p2 <- plot_branch_composition(x, style = "stack", by = "sample",
                                palette = palettes$composition) +
    ggplot2::labs(title = "nosZ branch partition")
  p3 <- try(plot_alpha_box(x, index = "shannon", branch = branch,
                           palette = palettes$alpha, facet = FALSE) +
              ggplot2::labs(title = paste0(branch, " Shannon diversity")),
            silent = TRUE)
  p4 <- try(plot_ordination(x, branch, method = "pcoa",
                            palette = palettes$beta, label = FALSE) +
              ggplot2::labs(title = paste0(branch, " PCoA (Bray-Curtis)")),
            silent = TRUE)
  if (inherits(p3, "try-error")) p3 <- patchwork::plot_spacer()
  if (inherits(p4, "try-error")) p4 <- patchwork::plot_spacer()
  # 显式经 patchwork 组合，避免依赖运行时已 attach 的 | 方法
  top <- patchwork::wrap_plots(p1, p2, ncol = 2)
  bottom <- patchwork::wrap_plots(p3, p4, ncol = 2)
  patchwork::wrap_plots(top, bottom, ncol = 1) +
    patchwork::plot_annotation(tag_levels = "A",
                               title = "NosZRef result overview",
                               theme = theme_nosz())
}

#' One-click visualization of a complete NosZRef run
#'
#' Generates the full figure collection (abundance, composition, alpha, beta,
#' differential abundance and network panels) for every available branch,
#' writes statistical tables, and builds a browsable HTML gallery.
#'
#' @param x A `noszref` object, or a path passed to [read_noszref()].
#' @param outdir Output directory.
#' @param formats Figure formats, any of `"pdf"`, `"png"`, `"tiff"`, `"svg"`.
#' @param dpi Raster resolution for png/tiff.
#' @param design Optional design passed to [read_noszref()] when `x` is a path.
#' @param levels Taxonomic levels used for composition figures.
#' @param diff_level Taxonomic level used for differential abundance.
#' @param permutations Permutations for PERMANOVA.
#' @param index Build an HTML index.
#' @param verbose Print progress.
#' @return Invisibly, a list with `files`, `failed` and `index`.
#' @export
noszref_visualize <- function(x, outdir = "nosZviz_output",
                              formats = c("pdf", "png"), dpi = 600,
                              design = NULL, env = NULL,
                              levels = c("phylum", "class",
                                                        "order", "family",
                                                        "genus", "species"),
                              diff_level = "genus", permutations = 999,
                              index = TRUE, verbose = TRUE) {
  if (is.character(x)) x <- read_noszref(x, design = design, verbose = verbose)
  stopifnot(inherits(x, "noszref"))
  env_ok <- !is.null(tryCatch(prepare_env(x, env), error = function(e) NULL))
  if (!is.null(env) && !env_ok) {
    nosz_info("  [warn] 环境表无法解析（缺 SampleID 或数值列），跳过环境关联图。")
  }
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  figdir <- file.path(outdir, "figures")
  tabdir <- file.path(outdir, "tables")
  dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)

  files <- character(0)
  failed <- list()
  gallery <- list()
  primary_b <- .primary_branch(x)

  run_one <- function(id, title, note, expr, width = 89, height = 59,
                      category = "General", code = NULL) {
    f <- file.path(figdir, id)
    res <- try(suppressWarnings(eval(expr)), silent = TRUE)
    if (inherits(res, "try-error") || is.null(res)) {
      failed[[id]] <<- as.character(res)
      if (verbose) nosz_info("  [skip] ", id, " : ",
                             gsub("\n", " ", as.character(res)))
      return(invisible(NULL))
    }
    ok <- try(nosz_save(res, f, width = width, height = height,
                        formats = formats, dpi = dpi), silent = TRUE)
    if (inherits(ok, "try-error")) {
      failed[[id]] <<- as.character(ok)
      if (verbose) nosz_info("  [fail] ", id, " : ", as.character(ok))
      return(invisible(NULL))
    }
    files <<- c(files, ok)
    # 自动推导展示代码：把调用反序列化，x 换成 res；过长则回落到模板
    if (is.null(code)) {
      code <- paste(deparse(expr, width.cutoff = 400), collapse = " ")
      code <- gsub("\\s+", " ", code)
      code <- gsub("\\bx\\b", "res", code)
      if (nchar(code) > 240) code <- ""
    }
    gallery[[length(gallery) + 1]] <<- list(
      id = id, title = title, note = note, category = category,
      code = code,
      png = paste0(id, ".png"), pdf = paste0(id, ".pdf")
    )
    if (verbose) nosz_info("  [ok]   ", id)
    invisible(NULL)
  }

  # ---------------- 全局图 ----------------
  run_one("01_branch_rpkm_bar", "Branch abundance (bar)",
          "各样本 nosZ-I/II/III 总 RPKM", quote(plot_branch_rpkm(x)),
          89, 62, "Abundance")
  run_one("02_branch_rpkm_summary", "Branch abundance (group summary)",
          "组均值 ± SD，散点为样本", quote(plot_branch_rpkm(x, style = "summary")),
          89, 62, "Abundance")
  run_one("03_branch_composition_stack", "Branch partition (stacked)",
          "三分支相对占比（按样本）",
          quote(plot_branch_composition(x, style = "stack")), 89, 62, "Abundance")
  run_one("04_branch_composition_donut", "Branch partition (donut)",
          "组均值三分支占比环形图",
          quote(plot_branch_composition(x, style = "donut")), 120, 55, "Abundance")
  run_one("05_reads_partition", "Read partition (Sankey)",
          "原始 reads 到三分支的流向",
          quote(plot_reads_partition(x)), 110, 65, "Abundance")
  run_one("06_branch_violin", "Branch abundance distribution",
          "分支丰度分布小提琴图", quote(plot_branch_violin(x)), 89, 62, "Abundance")
  run_one("07_branch_heatmap", "Branch abundance heatmap",
          "样本 × 分支 z-score 热图", quote(plot_branch_heatmap(x)), 100, 62, "Abundance")
  run_one("08_rank_abundance", "Rank-abundance curves",
          "OTU 秩-丰度曲线", quote(plot_rank_abundance(x)), 89, 62, "Abundance")
  run_one("09_overview_dashboard", "Overview dashboard",
          "四面板总览", quote(plot_overview_dashboard(x)), 183, 130, "Overview")
  run_one("10_tax_ternary", "Ternary plot of branch contribution",
          "属水平三分支贡献三角图", quote(plot_tax_ternary(x)), 110, 95, "Abundance")

  # ---------------- 支系生态学叙事（新增模块） ----------------
  run_one("11_clade_profile_phylum", "Clade-resolved composition (phylum)",
          "支系嵌套门水平组成：支系段描边 + 门渐变填充",
          quote(plot_clade_profile(x, level = "phylum", top_n = 6)),
          120, 80, "Clade")
  run_one("12_clade_profile_genus", "Clade-resolved composition (genus)",
          "支系嵌套属水平组成",
          quote(plot_clade_profile(x, level = "genus", top_n = 8)),
          120, 85, "Clade")
  run_one("13_clade_profile_group", "Clade composition by group",
          "组均值支系嵌套组成",
          quote(plot_clade_profile(x, level = "phylum", top_n = 6, by = "group")),
          100, 75, "Clade")
  if (all(c("nosZI", "nosZII") %in% x$branches)) {
    run_one("14_clade_ratio", "Clade II / Clade I ratio",
            "log2(nosZII/nosZI) 比值柱 + 组均值虚线",
            quote(plot_clade_ratio(x)), 110, 70, "Clade")
  }
  run_one("15_clade_scatter", "Pairwise clade abundance scatter",
            "支系间 RPKM 两两散点（Spearman）",
          quote(plot_clade_scatter(x)), 150, 70, "Clade")
  utils::write.csv(calc_clade_stats(x),
                   file.path(tabdir, "clade_stats.csv"), row.names = FALSE)

  # ---------------- 公众号风格扩展图型 ----------------
  run_one("16_chord_group_clade", "Group-clade chord diagram",
          "分组 <-> 支系丰度流向弦图（circlize）",
          quote(plot_chord(x, value = "relative")), 130, 120, "Clade")
  run_one("17_circular_bar", "Circular clade bar plot",
          "分组 x 支系环形柱状图",
          quote(plot_circular_bar(x)), 120, 105, "Clade")
  run_one("18_circular_heatmap", "Circular sample-taxon heatmap",
          "样本 x 类群 z-score 环形热图（按支系分扇区）",
          quote(plot_circular_heatmap(x)), 120, 120, "Composition")
  run_one("19_taxa_bump", "Taxa rank bump chart",
          "Top 类群丰度排名随分组变化",
          bquote(plot_taxa_bump(x, .(primary_b))), 130, 80, "Differential")

  # ---------------- 环境关联图（提供 env 时） ----------------
  if (env_ok) {
    run_one("20_envfit_pcoa", "PCoA with environmental fit",
            "PCoA + envfit 环境向量（R2/P 标注）",
            bquote(plot_envfit_ordination(x, .(primary_b), env = .(env))),
            110, 85, "Environment",
            code = sprintf('plot_envfit_ordination(res, "%s", env = env)', primary_b))
    run_one("21_alpha_env_corr", "Alpha x environment correlation",
            "alpha 指数与环境因子 Spearman 热图",
            bquote(plot_alpha_env_corr(x, .(primary_b), env = .(env))),
            100, 65, "Environment",
            code = sprintf('plot_alpha_env_corr(res, "%s", env = env)', primary_b))
    run_one("22_mantel_combo", "Mantel test combo",
            "环境 Spearman 热图 + 分支群落 Mantel 连线（公众号经典组合图）",
            bquote(plot_mantel_env(x, env = .(env))), 160, 95, "Environment",
            code = "plot_mantel_env(res, env = env)")
    run_one("23_env_corrplot", "Environment correlation heatmap",
            "环境因子相关热图（数值 + 显著性星号）",
            bquote(plot_env_corrplot(x, env = .(env))), 100, 95, "Environment",
            code = "plot_env_corrplot(res, env = env)")
    for (b in x$branches) {
      if (is.null(x$tax[[b]])) next
      run_one(sprintf("%s_rda", b), paste0(b, " RDA"),
              "约束排序 RDA（类群 x 环境因子）",
              bquote(plot_ordination_constrained(x, .(b), env = .(env),
                                                 method = "RDA")),
              110, 85, "Environment",
              code = sprintf('plot_ordination_constrained(res, "%s", env = env, method = "RDA")', b))
      run_one(sprintf("%s_cca", b), paste0(b, " CCA"),
              "对应分析 CCA（单峰模型）",
              bquote(plot_ordination_constrained(x, .(b), env = .(env),
                                                 method = "CCA")),
              110, 85, "Environment",
              code = sprintf('plot_ordination_constrained(res, "%s", env = env, method = "CCA")', b))
      run_one(sprintf("%s_taxa_env_network", b), paste0(b, " taxa-env network"),
              "类群-环境因子相关网络（Spearman）",
              bquote(plot_taxa_env_network(x, .(b), env = .(env))),
              130, 110, "Environment",
              code = sprintf('plot_taxa_env_network(res, "%s", env = env)', b))
    }
  }

  # ---------------- 分分支图 ----------------
  for (b in x$branches) {
    has_alpha <- !is.null(x$alpha[[b]])
    has_beta <- !is.null(x$beta[[b]])
    has_net <- !is.null(x$network[[b]])
    has_tax <- !is.null(x$tax[[b]])

    if (has_alpha) {
      run_one(sprintf("%s_alpha_panel", b), paste0(b, " alpha diversity panel"),
              "六种 alpha 指数箱线图",
              bquote(plot_alpha_panel(x, .(b))), 183, 90, "Alpha")
      run_one(sprintf("%s_alpha_violin", b), paste0(b, " alpha violin"),
              "Shannon 小提琴图",
              bquote(plot_alpha_violin(x, "shannon", .(b))), 89, 62, "Alpha")
      run_one(sprintf("%s_alpha_heatmap", b), paste0(b, " alpha heatmap"),
              "样本 × 指数 z-score 热图",
              bquote(plot_alpha_heatmap(x, .(b))), 110, 62, "Alpha")
      run_one(sprintf("%s_alpha_scatter", b), paste0(b, " alpha scatter"),
              "Shannon vs Richness 相关",
              bquote(plot_alpha_scatter(x, .(b))), 89, 70, "Alpha")
      run_one(sprintf("%s_rarefaction", b), paste0(b, " rarefaction"),
              "稀释曲线（组均值 ± SD）",
              bquote(plot_rarefaction(x, .(b))), 89, 62, "Alpha")
      run_one(sprintf("%s_alpha_ridge", b), paste0(b, " alpha ridgeline"),
              "Shannon 山脊图（组间密度分布）",
              bquote(plot_alpha_ridge(x, .(b))), 100, 70, "Alpha")
    }
    if (has_tax) {
      for (lv in intersect(c("phylum", "genus"), levels)) {
        run_one(sprintf("%s_tax_stack_%s", b, lv),
                paste0(b, " ", lv, " composition"),
                paste0(lv, " 水平组成堆叠柱"),
                bquote(plot_tax_stack(x, .(lv), .(b), top_n = 10)),
                89, 75, "Composition")
      }
      run_one(sprintf("%s_tax_heatmap_genus", b), paste0(b, " genus heatmap"),
              "属水平 Top30 热图",
              bquote(plot_tax_heatmap(x, .(b), "genus", top_n = 30)),
              110, 110, "Composition")
      run_one(sprintf("%s_tax_bubble_genus", b), paste0(b, " genus bubble"),
              "属水平气泡图",
              bquote(plot_tax_bubble(x, .(b), "genus", top_n = 20)),
              100, 95, "Composition")
      run_one(sprintf("%s_tax_area_phylum", b), paste0(b, " phylum alluvial"),
              "门水平河流图",
              bquote(plot_tax_area(x, .(b), "phylum", top_n = 8)),
              100, 70, "Composition")
      run_one(sprintf("%s_tax_sunburst", b), paste0(b, " sunburst"),
              "门-属双层环图",
              bquote(plot_tax_sunburst(x, .(b))), 130, 70, "Composition")
      run_one(sprintf("%s_top_taxa_genus", b), paste0(b, " top taxa by group"),
              "Top12 属分组条形图",
              bquote(plot_top_taxa(x, .(b), "genus", top_n = 12)),
              183, 110, "Composition")
      run_one(sprintf("%s_tax_corr_genus", b), paste0(b, " genus correlation"),
              "属水平 Spearman 相关",
              bquote(plot_tax_corr(x, .(b), "genus", top_n = 20)),
              120, 110, "Composition")
      # 支系生态学新增：类群树状图（丰度谱系 UPGMA）
      run_one(sprintf("%s_tax_dendrogram", b), paste0(b, " taxa dendrogram"),
              "类群丰度谱系 UPGMA 环形树（门着色、丰度定大小）",
              bquote(plot_tax_dendrogram(x, .(b), "genus", top_n = 40)),
              130, 115, "Clade")
    }
    # 组间共享 OTU（UpSet）
    if (!is.null(x$otutab[[b]]) && length(unique(x$design$Group)) >= 2) {
      run_one(sprintf("%s_otu_upset", b), paste0(b, " shared OTU (UpSet)"),
              "组间共享/特有 OTU 集合交叉",
              bquote(plot_group_upset(x, .(b), "otu")), 110, 80, "Differential")
    }
    # LEfSe 风格差异 + 指示种
    if (has_tax && length(unique(x$design$Group)) >= 2) {
      lf <- try(calc_lefse(x, b, diff_level, wilcox_p = 0.1, lda_cut = 0.8),
                silent = TRUE)
      if (!inherits(lf, "try-error") && nrow(lf) > 0) {
        utils::write.csv(lf, file.path(tabdir, paste0(b, "_lefse_", diff_level, ".csv")),
                         row.names = FALSE)
        run_one(sprintf("%s_lefse", b), paste0(b, " LEfSe-style"),
                "LEfSe 风格效应量条形图（KW + Wilcoxon）",
                bquote(plot_lefse(x, res = .(lf))),
                110, max(60, min(160, 20 + 6 * nrow(lf))), "Differential",
                code = sprintf('lf <- calc_lefse(res, "%s", "%s"); plot_lefse(res, res = lf)',
                               b, diff_level))
      }
      ind <- try(calc_indicator(x, b, diff_level, nperm = 499), silent = TRUE)
      if (!inherits(ind, "try-error") && nrow(ind) > 0) {
        utils::write.csv(ind,
                         file.path(tabdir, paste0(b, "_indicator_", diff_level, ".csv")),
                         row.names = FALSE)
        run_one(sprintf("%s_indicator", b), paste0(b, " indicator taxa"),
                "指示种 IndVal 棒棒糖图（A x B + 置换检验）",
                bquote(plot_indicator(x, res = .(ind))), 110, 90, "Differential",
                code = sprintf('ind <- calc_indicator(res, "%s", "%s"); plot_indicator(res, res = ind)',
                               b, diff_level))
      }
    }
    if (has_beta) {
      run_one(sprintf("%s_ordination_pcoa", b), paste0(b, " PCoA"),
              "Bray-Curtis PCoA + PERMANOVA",
              bquote(plot_ordination(x, .(b), "pcoa")), 100, 78, "Beta")
      run_one(sprintf("%s_ordination_nmds", b), paste0(b, " NMDS"),
              "NMDS 排序",
              bquote(plot_ordination(x, .(b), "nmds")), 100, 78, "Beta")
      run_one(sprintf("%s_ordination_pca", b), paste0(b, " PCA"),
              "Hellinger PCA",
              bquote(plot_ordination(x, .(b), "pca")), 100, 78, "Beta")
      run_one(sprintf("%s_beta_panel", b), paste0(b, " multi-distance PCoA"),
              "四种距离矩阵 PCoA 对比",
              bquote(plot_beta_panel(x, .(b))), 150, 130, "Beta")
      run_one(sprintf("%s_beta_heatmap", b), paste0(b, " distance heatmap"),
              "距离矩阵热图",
              bquote(plot_beta_heatmap(x, .(b))), 110, 95, "Beta")
      run_one(sprintf("%s_beta_dendrogram", b), paste0(b, " UPGMA tree"),
              "UPGMA 聚类树 + 分组",
              bquote(plot_beta_dendrogram(x, .(b))), 120, 80, "Beta")
      run_one(sprintf("%s_beta_dispersion", b), paste0(b, " beta dispersion"),
              "组内离散度（betadisper）",
              bquote(plot_beta_dispersion(x, .(b))), 89, 62, "Beta")
      run_one(sprintf("%s_sample_corr", b), paste0(b, " sample correlation"),
              "样本 Spearman 相关热图",
              bquote(plot_sample_corr(x, .(b))), 100, 90, "Beta")
    }
    if (has_alpha && has_beta && length(unique(x$design$Group)) >= 2) {
      run_one(sprintf("%s_diff_volcano", b), paste0(b, " differential volcano"),
              "差异属火山图",
              bquote(plot_diff_volcano(x, NULL, .(b), .(diff_level))),
              100, 80, "Differential")
      run_one(sprintf("%s_diff_ma", b), paste0(b, " MA plot"),
              "MA 图", bquote(plot_diff_ma(x, NULL, .(b), .(diff_level))),
              100, 75, "Differential")
      run_one(sprintf("%s_diff_forest", b), paste0(b, " effect size forest"),
              "效应量森林图",
              bquote(plot_diff_forest(x, NULL, .(b), .(diff_level), "t.test")),
              100, 85, "Differential")
      run_one(sprintf("%s_diff_stamp", b), paste0(b, " STAMP-style"),
              "STAMP 扩展误差条图",
              bquote(plot_diff_stamp(x, NULL, .(b), .(diff_level), "t.test")),
              160, 85, "Differential")
      run_one(sprintf("%s_diff_heatmap", b), paste0(b, " differential heatmap"),
              "差异属热图（含显著性星号）",
              bquote(plot_diff_heatmap(x, NULL, .(b), .(diff_level))),
              110, 95, "Differential")
    }
    if (has_net) {
      for (g in unique(x$network[[b]]$nodes$Group)) {
        run_one(sprintf("%s_network_%s", b, g),
                paste0(b, " network - ", g),
                paste0(g, " 组共现网络（模块着色）"),
                bquote(plot_network(x, .(b), .(g))), 130, 110, "Network")
      }
      run_one(sprintf("%s_network_zipi", b), paste0(b, " Zi-Pi"),
              "节点拓扑角色 Zi-Pi 图",
              bquote(plot_network_zipi(x, .(b))), 120, 85, "Network")
      run_one(sprintf("%s_network_properties", b), paste0(b, " network properties"),
              "网络拓扑参数对比",
              bquote(plot_network_properties(x, .(b))), 183, 90, "Network")
      run_one(sprintf("%s_network_degree", b), paste0(b, " degree distribution"),
              "度分布（log-log）",
              bquote(plot_network_degree(x, .(b))), 89, 70, "Network")
      run_one(sprintf("%s_network_modules", b), paste0(b, " network modules"),
              "模块规模与连通性",
              bquote(plot_network_modules(x, .(b))), 140, 85, "Network")
      run_one(sprintf("%s_network_upset", b), paste0(b, " shared nodes"),
              "组间共享节点 UpSet 图",
              bquote(plot_network_upset(x, .(b))), 120, 85, "Network")
    }
  }

  # ---------------- 统计表 ----------------
  utils::write.csv(noszref_summary(x), file.path(tabdir, "nosz_summary.csv"),
                   row.names = FALSE)
  if (!is.null(x$design)) {
    utils::write.csv(x$design, file.path(tabdir, "design.csv"), row.names = FALSE)
  }
  for (b in x$branches) {
    if (!is.null(x$beta[[b]]) && length(unique(x$design$Group)) >= 2) {
      dists <- names(x$beta[[b]])
      st <- do.call(rbind, lapply(dists, function(dd) {
        try(calc_permanova(x, b, dd, permutations), silent = TRUE)
      }))
      st <- st[!vapply(st, inherits, logical(1), "try-error")]
      if (length(st) > 0) {
        utils::write.csv(do.call(rbind, st),
                         file.path(tabdir, paste0(b, "_beta_statistics.csv")),
                         row.names = FALSE)
      }
    }
    if (!is.null(x$network[[b]]) && !is.null(x$network[[b]]$stats)) {
      utils::write.csv(x$network[[b]]$stats,
                       file.path(tabdir, paste0(b, "_network_statistics.csv")),
                       row.names = FALSE)
    }
    if (!is.null(x$tax[[b]]) && length(unique(x$design$Group)) >= 2) {
      dd <- try(calc_diff_taxa(x, b, diff_level), silent = TRUE)
      if (!inherits(dd, "try-error")) {
        utils::write.csv(dd, file.path(tabdir, paste0(b, "_diff_", diff_level, ".csv")),
                         row.names = FALSE)
      }
    }
  }

  idx <- NA_character_
  if (index) {
    idx <- file.path(outdir, "index.html")
    write_gallery(gallery, idx, x)
  }
  if (verbose) {
    nosz_info("完成：", length(files), " 个文件；失败/跳过 ",
              length(failed), " 项。输出目录：", normalizePath(outdir))
  }
  invisible(list(files = files, failed = failed, index = idx,
                 gallery = gallery))
}

# 生成 HTML 图集
write_gallery <- function(gallery, path, x) {
  cats <- unique(vapply(gallery, `[[`, character(1), "category"))
  css <- paste(
    "body{font-family:-apple-system,'Segoe UI','Microsoft YaHei',sans-serif;",
    "margin:0;padding:22px 28px;background:#f6f7f9;color:#1c2530;}",
    "h1{font-size:20px;margin:0 0 4px;}h2{font-size:14px;margin:26px 0 10px;",
    "padding-bottom:4px;border-bottom:1px solid #d7dde5;color:#33475f;}",
    ".meta{font-size:12px;color:#5b6b7c;margin-bottom:18px;}",
    ".grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));",
    "gap:14px;}",".card{background:#fff;border:1px solid #e3e8ee;border-radius:6px;",
    "padding:8px;box-shadow:0 1px 2px rgba(20,40,70,.05);}",
    ".card img{width:100%;height:auto;display:block;border-radius:3px;",
    "background:#fff;}",
    ".cap{font-size:11.5px;color:#3d4d5e;margin-top:6px;line-height:1.45;}",
    ".cap b{color:#12202e;font-size:12px;}.tag{display:inline-block;font-size:10px;",
    "color:#456;background:#eef2f7;border-radius:3px;padding:1px 5px;margin-left:4px;}",
    "a{color:#2b6cb0;text-decoration:none;font-size:10.5px;}",
    "pre.code{background:#f2f5f8;border:1px solid #e0e6ec;border-radius:4px;",
    "padding:5px 7px;margin:5px 0 0;font-size:10px;color:#1a3a5c;",
    "overflow-x:auto;white-space:pre-wrap;word-break:break-all;}",
    "button.cp{float:right;font-size:9.5px;border:1px solid #cdd6e0;",
    "background:#fff;border-radius:3px;padding:1px 6px;cursor:pointer;",
    "color:#456;}",
    "button.cp:hover{background:#eef2f7;}"
  )
  head <- c(
    "<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'>",
    "<title>NosZRef visualization gallery</title>",
    paste0("<style>", paste(css, collapse = ""), "</style>"),
    paste0("<script>function cp(t){var e=t.parentElement.querySelector('pre');",
           "navigator.clipboard.writeText(e.innerText);t.textContent='已复制';",
           "setTimeout(function(){t.textContent='复制代码'},1200)}</script></head><body>"),
    "<h1>NosZRef 可视化图库 / Visualization gallery</h1>",
    paste0("<div class='meta'>由 nosZviz ", utils::packageVersion("nosZviz"),
           " 生成 &nbsp;|&nbsp; 分支：", paste(x$branches, collapse = " / "),
           " &nbsp;|&nbsp; 样本：", nrow(x$design),
           " &nbsp;|&nbsp; 分组：", paste(unique(x$design$Group), collapse = " / "),
           " &nbsp;|&nbsp; 共 ", length(gallery), " 张图",
           "<br>用法：找到想要的图 → 复制卡片下方的一行代码 → ",
           "把 <code>res</code> 换成你的 <code>read_noszref()</code> 结果即可",
           "（环境类图还需 <code>env</code> 环境因子表，含 SampleID 列）。",
           "演示图所用环境数据为模拟值。</div>")
  )
  body <- unlist(lapply(cats, function(ct) {
    items <- Filter(function(g) identical(g$category, ct), gallery)
    c(paste0("<h2>", ct, " (", length(items), ")</h2>"),
      c("<div class='grid'>", unlist(lapply(items, function(g) {
        code_block <- ""
        if (!is.null(g$code) && nzchar(g$code)) {
          code_block <- paste0(
            "<button class='cp' onclick='cp(this)'>复制代码</button>",
            "<pre class='code'>", gsub("<", "&lt;", g$code, fixed = TRUE),
            "</pre>")
        }
        c("<div class='card'>",
          paste0("<a href='figures/", g$pdf, "' target='_blank'>",
                 "<img src='figures/", g$png, "' alt='", g$id, "'></a>"),
          paste0("<div class='cap'><b>", g$title, "</b><span class='tag'>",
                 g$id, "</span><br>", g$note, "<br><a href='figures/",
                 g$pdf, "' target='_blank'>PDF</a>",
                 code_block, "</div>"),
          "</div>")
      })), "</div>"))
  }))
  html <- c(head, body, "</body></html>")
  # 用二进制连接 + 只纠正编码标记的方式写出，避免非 UTF-8 locale 下
  # enc2utf8() 把已正确的 UTF-8 字节二次编码（中文变乱码）
  nosz_write_utf8(html, path)
  invisible(path)
}
