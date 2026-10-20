#' nosZviz: publication-ready visualization for NosZRef results
#'
#' The package reads the complete output directory of the NosZRef pipeline
#' and provides a coherent grammar of graphics for nosZ-I / nosZ-II / nosZ-III
#' profiling: gene abundance partitioning, taxonomic composition, alpha and
#' beta diversity, differential abundance and co-occurrence networks.
#'
#' @keywords internal
#' @importFrom dplyr %>% .data
#' @importFrom magrittr %>%
#' @importFrom rlang .data
"_PACKAGE"

utils::globalVariables(c(".", "Abundance", "Branch", "Combo", "Degree",
                         "Degree2", "Freq", "Group", "Hub", "In", "Index",
                         "Mean", "Metric", "Module", "Node", "P", "Pi",
                         "Rank", "Role", "SampleID", "Set", "Taxon", "Unit",
                         "Value", "Zi", "assigned", "ci_high", "ci_low",
                         "f", "frac", "key", "label", "leaf", "level1",
                         "level2", "m", "mean1", "mean2", "name", "node",
                         "p", "padj", "prevalence", "px", "py", "r", "s",
                         "sd", "sig", "stratum", "v", "value", "x", "y",
                         "ymax", "ymin"))

# ---------------------------------------------------------------------------
# 小工具：默认值、报错、静默读取
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# 统一风格的报错
nosz_abort <- function(...) {
  stop(paste0(...), call. = FALSE)
}

# 统一风格的信息输出
nosz_info <- function(...) {
  message(paste0(...))
}

# 判断某个分支的某个结果是否存在，不存在时给出清晰提示
nosz_need <- function(obj, what, branch = NULL) {
  x <- obj[[what]]
  if (is.null(branch)) {
    if (is.null(x)) {
      nosz_abort("NosZRef 结果中缺少 `", what, "`（可能该步骤未完成）。")
    }
    return(x)
  }
  if (is.null(x) || is.null(x[[branch]])) {
    nosz_abort("NosZRef 结果中缺少分支 `", branch, "` 的 `", what,
               "`。可用分支：", paste(obj$branches, collapse = ", "))
  }
  x[[branch]]
}

# 读取分隔符文本时容忍空文件与 BOM
read_delim_safe <- function(file, delim = "\t", ...) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  if (file.info(file)$size == 0) return(NULL)
  out <- try(
    readr::read_delim(file, delim = delim, show_col_types = FALSE,
                      progress = FALSE, ...),
    silent = TRUE
  )
  if (inherits(out, "try-error")) return(NULL)
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

# 数据框 -> 数值矩阵（第一列作行名）
df_to_matrix <- function(df, id_col = 1) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  ids <- as.character(df[[id_col]])
  m <- as.matrix(df[, -id_col, drop = FALSE])
  m <- apply(m, 2, function(v) suppressWarnings(as.numeric(v)))
  if (is.null(dim(m))) m <- matrix(m, nrow = length(ids))
  rownames(m) <- ids
  colnames(m) <- colnames(df)[-id_col]
  m[is.na(m)] <- 0
  m
}

# 长表 -> 宽矩阵（taxon x sample）
long_to_matrix <- function(df, tax_col, sample_col, value_col) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df[, c(tax_col, sample_col, value_col)]
  colnames(df) <- c("taxon", "sample", "value")
  df$value <- suppressWarnings(as.numeric(df$value))
  df$value[is.na(df$value)] <- 0
  out <- stats::xtabs(value ~ taxon + sample, data = df)
  matrix(as.numeric(out), nrow = nrow(out),
         dimnames = list(rownames(out), colnames(out)))
}

# 矩阵 -> 长表，并接上分组信息
matrix_to_long <- function(mat, design, group_col = "Group") {
  if (is.null(mat) || nrow(mat) == 0) return(NULL)
  df <- data.frame(
    Taxon = rep(rownames(mat), times = ncol(mat)),
    SampleID = rep(colnames(mat), each = nrow(mat)),
    Abundance = as.numeric(mat),
    stringsAsFactors = FALSE
  )
  if (!is.null(design)) df <- merge(df, design, by = "SampleID", all.x = TRUE)
  df
}

# 行标准化（z-score），用于热图
zscore_rows <- function(mat) {
  if (is.null(mat) || nrow(mat) == 0) return(mat)
  z <- t(apply(mat, 1, function(v) {
    s <- stats::sd(v)
    if (is.na(s) || s == 0) return(rep(0, length(v)))
    (v - mean(v)) / s
  }))
  dimnames(z) <- dimnames(mat)
  z
}

# 相对丰度（列内归一，0-1）
rel_abund <- function(mat) {
  cs <- colSums(mat, na.rm = TRUE)
  cs[cs == 0] <- 1
  sweep(mat, 2, cs, "/")
}

# 合并低丰度分类单元为 "Other"
collapse_taxa <- function(mat, top_n = 10, by = "mean") {
  if (is.null(mat) || nrow(mat) == 0) return(mat)
  if (nrow(mat) <= top_n) return(mat)
  score <- if (identical(by, "mean")) rowMeans(mat) else rowSums(mat)
  ord <- order(score, decreasing = TRUE)
  keep <- ord[seq_len(top_n)]
  rest <- setdiff(seq_len(nrow(mat)), keep)
  out <- rbind(
    mat[keep, , drop = FALSE],
    Other = colSums(mat[rest, , drop = FALSE])
  )
  out[order(rowSums(out), decreasing = TRUE), , drop = FALSE]
}

# 解析 "d_Bacteria;p_Pseudomonadota;..." 形式的分类字符串
tax_parse <- function(x, prefix = c("d", "p", "c", "o", "f", "g", "s")) {
  x <- as.character(x)
  n <- length(x)
  parts <- strsplit(x, ";", fixed = TRUE)
  get_rank <- function(pat) {
    vapply(parts, function(p) {
      hit <- grep(pat, p, value = TRUE)
      if (length(hit) == 0) return(NA_character_)
      sub(paste0("^", pat, "_+"), "", hit[1])
    }, character(1))
  }
  out <- data.frame(
    Kingdom = get_rank("d|k"), Phylum = get_rank("p"), Class = get_rank("c"),
    Order = get_rank("o"), Family = get_rank("f"), Genus = get_rank("g"),
    Species = get_rank("s"), stringsAsFactors = FALSE
  )
  # 有些版本使用 k__ 前缀，统一到 Kingdom
  out
}

# 组间配色：优先使用手工调色板，超出则自动扩展
group_palette <- function(groups, palette = "nosz") {
  groups <- unique(as.character(groups))
  n <- length(groups)
  pal <- nosz_pal(n, palette)
  stats::setNames(pal, groups)
}

# 取第一个非空值
first_nonempty <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0) NA_character_ else x[1]
}

# 安全的 p 值格式化
fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = digits)))
}

# 星号标记
p_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# ---------------------------------------------------------------------
# UTF-8 文本写出工具（解决中文乱码）
#
# 背景：Windows 下若 R 的 locale 不是 UTF-8（例如 locale = C，或
# "Setting LC_CTYPE=C.UTF-8 failed"），R 会把 UTF-8 源码文本按 latin1
# 解析：**字节仍是合法 UTF-8，但被错标成 "latin1"**。
# 此时调用 enc2utf8() 会把每个字节再编码一次，中文就变成
# "Ã¥Â¯Ã¨Â§â€ ..." 这样的乱码（index.html 图库曾出现此问题）。
#
# 对策：先判断字节本身是否为合法 UTF-8——是则只纠正编码标记、不动字节；
# 否则（确实是本地编码如 GBK）再做正常转换。最后用二进制连接写出，
# 彻底绕开 locale 带来的隐式转换。
# ---------------------------------------------------------------------

# 把字符向量规范成「标记为 UTF-8」的形态
nosz_utf8 <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return(x)
  if (all(validUTF8(x))) {
    # 字节已是合法 UTF-8，只是标记错了：仅纠正标记
    Encoding(x) <- "UTF-8"
  } else {
    # 真正的本地编码（GBK/GB18030 等）走正常转换
    x <- enc2utf8(x)
    ok <- validUTF8(x)
    if (any(ok)) Encoding(x[ok]) <- "UTF-8"
  }
  x
}

# 以 UTF-8 字节写出文本文件（不做任何 locale 转换）
nosz_write_utf8 <- function(x, path, ...) {
  con <- file(path, open = "wb")          # 二进制连接：字节原样落盘
  on.exit(close(con), add = TRUE)
  writeLines(nosz_utf8(x), con, useBytes = TRUE, ...)
  invisible(path)
}


# 让每个分支的样本顺序与 design 一致
order_samples <- function(samples, design) {
  if (is.null(design)) return(samples)
  d <- design$SampleID[design$SampleID %in% samples]
  c(d, setdiff(samples, d))
}
