# nosZviz 安装教程（中文）

> 对应版本 v0.2.0。完整图文版见仓库 `workspace/nosZviz_安装教程.html`；
> 包内可用 `system.file("INSTALL_zh.md", package = "nosZviz")` 定位本文件。

## 1. 环境要求

| 项目 | 要求 |
|---|---|
| R 版本 | ≥ 4.1.0（建议 4.3+） |
| Rtools | **不需要**（纯 R 包，无 C/C++ 源码） |
| 网络 | 仅首次装依赖需要联网 |
| 系统 | Windows / macOS / Linux |

```r
R.version.string      # 确认版本
```

## 2. 安装依赖（只需一次）

```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

install.packages(c(
  "ggplot2", "dplyr", "magrittr", "rlang", "tidyr", "tibble", "readr",
  "purrr", "stringr", "forcats", "scales",
  "vegan", "igraph", "ggraph", "patchwork",
  "RColorBrewer", "ggrepel", "pheatmap", "ggpubr", "ggalluvial",
  "UpSetR", "ggridges", "circlize", "viridisLite", "gridExtra"
))
```

依赖用途：`vegan` 负责 beta/环境层全部统计；`igraph`+`ggraph` 负责网络；
`circlize` 负责弦图/环形图；`UpSetR` 负责 UpSet；`ggridges` 负责山脊图；
`ggalluvial` 负责桑基/河流图。缺哪个只会影响对应图型。

## 3. 安装 nosZviz

```r
# A. 从 tar.gz（推荐）
install.packages("D:/R/nosZviz_0.2.0.tar.gz", repos = NULL, type = "source")

# B. 从源码目录（含 DESCRIPTION 的 nosZviz/ 文件夹）
install.packages("D:/R/nosZviz", repos = NULL, type = "source")
devtools::install_local("D:/R/nosZviz", upgrade = "never")   # 或用 devtools

# C. 从 GitHub（仓库上线后）
remotes::install_github("ZhangBaoshan668/nosZviz", build_vignettes = TRUE)

# D. RStudio：Tools → Install Packages → Package Archive File → 选中 tar.gz
```

`repos = NULL` 与 `type = "source"` 两个参数都不能省。

## 4. 验证

```r
library(nosZviz)
packageVersion("nosZviz")     # 应为 0.2.0
find.package("nosZviz")       # 安装位置

fns <- c("read_noszref", "noszref_visualize", "plot_clade_profile",
         "plot_mantel_env", "plot_chord", "plot_ordination_constrained",
         "plot_lefse", "plot_alpha_ridge", "nosz_save", "theme_nosz")
sapply(fns, exists, where = asNamespace("nosZviz"))   # 全部 TRUE 即正常
```

## 5. 跑起来

```r
res <- read_noszref("D:/project/NosZRef/Example_data")   # 含 merge/ 的那层目录
plot_branch_rpkm(res)                                    # 单张图

noszref_visualize(res, outdir = "nosZviz_output",        # 一键全套
                  formats = c("pdf", "png"), dpi = 600)  # → 打开 index.html

# 有环境因子表时（CSV，含 SampleID 列 + 数值列）
noszref_visualize(res, outdir = "out_env", env = "env.csv")
```

命令行：`Rscript run_noszref_visual.R <输入目录> <输出目录> [metadata.txt] [env.csv]`

## 6. 常见报错

| 现象 | 解决 |
|---|---|
| `there is no package called 'xxx'` | 依赖没装全，回第 2 步或单独装该包 |
| `cannot open URL ...tar.gz` | 漏 `repos = NULL`，或路径不对（用 `/` 或 `\\`） |
| `had non-zero exit status` | 看上方日志，多为缺依赖；目录用方式 B、tar.gz 用方式 A |
| `not available for this version of R` | 换镜像 `https://cloud.r-project.org`，仍不行则升级 R |
| `Permission denied` | 默认库无写权限：`install.packages(..., lib="D:/Rlibs")` + `.libPaths("D:/Rlibs")` |
| 装完 `library()` 找不到 | 装到了非默认库：`library(nosZviz, lib.loc="D:/Rlibs")` |
| 图里中文变方框 | PDF 字体限制；包内图文字为英文。自定义中文标题建议导出 PNG 或配 `showtext` |
| 图库 HTML 中文乱码 | 旧版 bug，0.2.0 已修复；确认版本并重装 |
| `could not find function "plot_clade_profile"` | 装的是 0.1.0：`remove.packages("nosZviz")` 后重装 0.2.0 |
| 某分支图全部跳过 | 该分支 OTU 太少（如示例 nosZIII 仅 1 个），无法计算，属预期 |

## 7. 升级 / 卸载 / 离线

```r
install.packages("D:/R/nosZviz_0.2.0.tar.gz", repos = NULL, type = "source")  # 覆盖升级
remove.packages("nosZviz")                                                    # 卸载
remove.packages("nosZviz", lib = "D:/Rlibs")                                  # 指定库卸载
```

离线安装：联网机器上 `download.packages(pkgs, destdir="D:/offline", type="win.binary")`，
拷贝到目标机器后 `install.packages(list.files("D:/offline", pattern="\\.zip$", full.names=TRUE), repos=NULL, type="win.binary")`。
