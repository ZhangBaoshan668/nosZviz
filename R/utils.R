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
#' @importFrom utils head
#' @importFrom stats na.omit
#' @importFrom grDevices colorRampPalette
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
# 
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# 
nosz_abort <- function(...) {
  stop(paste0(...), call. = FALSE)
}

# 
nosz_info <- function(...) {
  message(paste0(...))
}

# 
nosz_need <- function(obj, what, branch = NULL) {
  x <- obj[[what]]
  if (is.null(branch)) {
    if (is.null(x)) {
      nosz_abort("NosZRef  `", what, "`")
    }
    return(x)
  }
  if (is.null(x) || is.null(x[[branch]])) {
    nosz_abort("NosZRef  `", branch, "`  `", what,
               "`", paste(obj$branches, collapse = ", "))
  }
  x[[branch]]
}

#  BOM
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

#  -> 
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

#  -> taxon x sample
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

#  -> 
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

# z-score
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

# 0-1
rel_abund <- function(mat) {
  cs <- colSums(mat, na.rm = TRUE)
  cs[cs == 0] <- 1
  sweep(mat, 2, cs, "/")
}

#  "Other"
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

#  "d_Bacteria;p_Pseudomonadota;..." 
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
  #  k__  Kingdom
  out
}

# 
group_palette <- function(groups, palette = "nosz") {
  groups <- unique(as.character(groups))
  n <- length(groups)
  pal <- nosz_pal(n, palette)
  stats::setNames(pal, groups)
}

# 
first_nonempty <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0) NA_character_ else x[1]
}

#  p 
fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = digits)))
}

# 
p_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# ---------------------------------------------------------------------
# UTF-8 
#
# Windows  R  locale  UTF-8 locale = C
# "Setting LC_CTYPE=C.UTF-8 failed"R  UTF-8  latin1
# ** UTF-8 "latin1"**
#  enc2utf8() 
# " ..." index.html 
#
#  UTF-8
#  GBK
#  locale 
# ---------------------------------------------------------------------

#  UTF-8
nosz_utf8 <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return(x)
  if (all(validUTF8(x))) {
    #  UTF-8
    Encoding(x) <- "UTF-8"
  } else {
    # GBK/GB18030 
    x <- enc2utf8(x)
    ok <- validUTF8(x)
    if (any(ok)) Encoding(x[ok]) <- "UTF-8"
  }
  x
}

#  UTF-8  locale 
nosz_write_utf8 <- function(x, path, ...) {
  con <- file(path, open = "wb")          # 
  on.exit(close(con), add = TRUE)
  writeLines(nosz_utf8(x), con, useBytes = TRUE, ...)
  invisible(path)
}


#  design 
order_samples <- function(samples, design) {
  if (is.null(design)) return(samples)
  d <- design$SampleID[design$SampleID %in% samples]
  c(d, setdiff(samples, d))
}
