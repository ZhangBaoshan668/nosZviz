#!/usr/bin/env Rscript
# =====================================================================
# nosZviz 一键可视化示例脚本
#
# 用法：
#   Rscript run_noszref_visual.R <NosZRef输出目录> <输出目录> [metadata.txt] [env.csv]
#
# 示例（NosZRef 自带示例数据）：
#   Rscript run_noszref_visual.R Example_data ./nosZviz_output
#   Rscript run_noszref_visual.R Example_data ./nosZviz_output "" env.csv
#
# 说明：
#   1) 自动解析 NosZRef 输出目录中的全部结果（merge/nosZI|II|III）；
#   2) 生成 60+ 张发表级图表（PDF + PNG）与统计表；
#   3) 生成 index.html 图库索引（每张图附一行可复制代码），浏览器打开即可浏览；
#   4) 提供环境因子表 env.csv 时（含 SampleID 列 + 数值列），
#      额外生成 RDA/CCA/db-RDA、envfit、Mantel 组合图、
#      alpha×环境相关热图、类群-环境相关网络等环境关联图。
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript run_noszref_visual.R <NosZRef输出目录> <输出目录> [metadata.txt] [env.csv]\n")
  quit(status = 1)
}

indir   <- args[1]                       # NosZRef 输出目录（含 merge/）
outdir  <- args[2]                       # 图表输出目录
design  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL  # 可选 metadata 文件
envf    <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NULL  # 可选环境因子表

library(nosZviz)

## 1. 读取数据 ---------------------------------------------------------
res <- read_noszref(indir, design = design, verbose = TRUE)
print(res)                               # 查看各分支结果可用性
print(noszref_summary(res))              # 概览统计表

## 2. 单图示例（可按需挑选） -------------------------------------------
# 分支丰度与占比
# print(plot_branch_rpkm(res))                       # 三分支 RPKM
# print(plot_branch_composition(res, style = "donut"))  # 三分支占比环形图
# 支系生态学（v0.2.0）
# print(plot_clade_profile(res, rank = "phylum"))    # 支系嵌套分类组成
# print(plot_clade_ratio(res))                       # Clade II/I log2 比值
# 环境因子关联（v0.2.0，需 env）
# print(plot_ordination_constrained(res, env, method = "RDA"))
# print(plot_mantel_env(res, env))                   # Mantel 组合图

## 3. 一键全套输出 ------------------------------------------------------
invisible(noszref_visualize(
  res,
  outdir      = outdir,          # 输出目录
  env         = envf,            # 环境因子表（NULL 则跳过环境层图表）
  formats     = c("pdf", "png"), # 每张图同时导出 PDF 与 PNG
  dpi         = 600,             # 位图分辨率
  levels      = c("phylum", "class", "order", "family", "genus", "species"),
  diff_level  = "genus",         # 差异分析使用的分类层级
  permutations = 999,            # PERMANOVA 置换次数
  index       = TRUE             # 生成 index.html 图库（带复制代码按钮）
))

cat("\n完成！请在浏览器中打开: ", file.path(outdir, "index.html"), "\n")
