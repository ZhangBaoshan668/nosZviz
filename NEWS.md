# nosZviz 0.2.0

## 新增：支系生态学叙事层（Clade）

- `plot_clade_profile()`：支系嵌套分类学组成图（支系段描边 + 支系内渐变填充），回答"Clade II 的优势由哪个类群驱动"
- `plot_clade_ratio()`：Clade II / Clade I log2 比值柱状图
- `plot_clade_scatter()`：支系两两 RPKM 散点 + Spearman 相关
- `calc_clade_stats()`：每样本支系丰度/占比统计表
- `plot_chord()`：分组 ↔ 支系丰度流向弦图（返回 `nosz_chord` S3 对象，`print()` 出图）
- `plot_circular_bar()` / `plot_circular_heatmap()`：环形柱状图 / 环形热图
- `plot_tax_dendrogram()`：类群丰度谱系 UPGMA 环形树（非系统发育树）

## 新增：环境因子关联层（Environment，需 `env=`）

- `prepare_env()`：环境表校验与样本对齐
- `plot_ordination_constrained()`：RDA / CCA / db-RDA（约束轴解释率 + 置换检验 P）
- `plot_envfit_ordination()`：PCoA + envfit 环境向量
- `plot_mantel_env()`：环境 Spearman 热图 + 支系群落 Mantel 连线组合图
- `plot_alpha_env_corr()`：alpha 指数 × 环境因子相关热图
- `plot_taxa_env_network()`：类群-环境因子相关网络
- `plot_env_corrplot()`：环境因子相关热图（原生 ggplot 实现，不再依赖 corrplot）

## 新增：差异进阶与其他

- `calc_lefse()` + `plot_lefse()`：LEfSe 风格效应量条形图
- `calc_indicator()` + `plot_indicator()`：指示种 IndVal（两组时完全枚举置换检验）
- `plot_group_upset()`：组间共享 / 特有 OTU UpSet
- `plot_taxa_bump()`：类群丰度排名 bump 图
- `plot_alpha_ridge()`：alpha 指数组间密度山脊图

## 改进

- `noszref_visualize()` 新增 `env=` 参数；图库新增 9 个分类块（v0.2.0 共 113 张图）
- 新增 `docs/` 图库站点（GitHub Pages）：完整 113 张图按 9 类分区、支持关键词搜索、
  每张图带「复制代码」，README 仅保留 12 张代表图（覆盖全部 9 个分析类型）
- HTML 图库（`index.html`）每张卡片附**一行可复制代码**与「复制代码」按钮
- 新增 `vignettes/gallery.Rmd` 图库文档与 `inst/INSTALL_zh.md` 安装教程
- `nosz_save()` 改用 `tryCatch` + `finally`，单个图报错不再污染后续图形的设备状态
- 修复 Windows 非 UTF-8 locale 下 HTML 图库中文乱码（新增 `nosz_write_utf8()`）
- `inst/scripts/run_noszref_visual.R` 支持第 4 个参数 `env.csv`

## 依赖变化

- 新增 Imports：`UpSetR`、`ggridges`、`circlize`
- 移除 Imports：`corrplot`

---

# nosZviz 0.1.0

首个版本。

- `read_noszref()`：解析 NosZRef 输出目录（分支 RPKM / reads、OTU 表、抽平表、7 级分类、alpha、6 种 beta 距离、BLAST 长表、网络拓扑）
- 覆盖丰度分配、分类组成、alpha/beta 多样性、差异丰度、共现网络五大类共 47 个绘图与统计函数
- `noszref_visualize()`：一键生成全套图表 + 统计表 + HTML 图集
- `theme_nosz()` / `nosz_pal()` / `nosz_save()`：统一主题与 mm 尺寸多格式导出
