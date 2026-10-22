# ---------------------------------------------------------------------------
# clade  NosZRef 
# Clade II/I  nosZ 
#  Jones et al. 2013 NatureHallin 
# ---------------------------------------------------------------------------

# 
clade_cols <- function(branches) {
  base <- c(nosZI = "#1F78B4", nosZII = "#FF7F00", nosZIII = "#33A02C")
  cols <- base[branches]
  if (length(cols) == 0) cols <- stats::setNames("#1F78B4", branches)
  cols[is.na(cols)] <- "grey55"
  cols
}

#  top_n  + 
clade_long <- function(x, level, top_n, by = c("sample", "group")) {
  by <- match.arg(by)
  out <- list()
  for (b in x$branches) {
    m <- try(nosz_tax_matrix(x, b, level, value = "rpkm"), silent = TRUE)
    if (inherits(m, "try-error") || is.null(m) || nrow(m) == 0) next
    m <- rel_abund(m)
    m <- collapse_taxa(m, top_n = top_n)
    if (by == "group") {
      gs <- unique(x$design$Group[x$design$SampleID %in% colnames(m)])
      m <- sapply(gs, function(g) {
        ss <- intersect(x$design$SampleID[x$design$Group == g], colnames(m))
        if (length(ss) == 0) return(rep(0, nrow(m)))
        rowMeans(m[, ss, drop = FALSE])
      })
      dimnames(m) <- list(rownames(m), gs)
    }
    df <- data.frame(Branch = b, Taxon = rep(rownames(m), times = ncol(m)),
                     ID = rep(colnames(m), each = nrow(m)),
                     Abundance = as.numeric(m), stringsAsFactors = FALSE)
    out[[length(out) + 1]] <- df
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

#' Clade-resolved nested taxonomic composition
#'
#' **Clade II
#' ** 
#' nosZ-I / II / III 
#' 
#'
#' @param x A `noszref` object.
#' @param level Taxonomic level nested within clades.
#' @param top_n Top taxa shown per clade, rest merged as "Other".
#' @param by `"sample"` `"group"` 
#' @param base_size Passed to [theme_nosz()].
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_clade_profile(res, level = "phylum", top_n = 6)
#' }
plot_clade_profile <- function(x, level = "phylum", top_n = 6,
                               by = c("sample", "group"), base_size = 8) {
  stopifnot(inherits(x, "noszref"))
  by <- match.arg(by)
  d <- clade_long(x, level, top_n, by = by)
  if (is.null(d)) nosz_abort(" ", level, " ")
  d$Branch <- factor(d$Branch, levels = x$branches)
  d$ID <- factor(d$ID, levels = order_samples(unique(d$ID), x$design))

  # 
  ord_tax <- names(sort(tapply(d$Abundance, paste(d$Branch, d$Taxon), sum),
                        decreasing = TRUE))
  d$TaxonB <- paste(d$Branch, d$Taxon)
  d$TaxonB <- factor(d$TaxonB, levels = ord_tax)

  #  -> 
  cols <- clade_cols(x$branches)
  fills <- unlist(lapply(x$branches, function(b) {
    tx <- unique(as.character(d$Taxon[d$Branch == b]))
    tx <- tx[order(match(paste(b, tx), ord_tax))]
    stats::setNames(grDevices::colorRampPalette(c(cols[[b]], "#F2F2F2"))(length(tx)),
                    paste(b, tx))
  }))
  lab_tax <- sub("^(nosZ[I]*) ", "", names(fills))
  names(lab_tax) <- names(fills)

  p <- ggplot2::ggplot(d, ggplot2::aes(ID, Abundance)) +
    ggplot2::geom_col(ggplot2::aes(fill = TaxonB, colour = Branch),
                      position = "stack", width = 0.72, linewidth = 0.55) +
    ggplot2::scale_fill_manual(values = fills, labels = lab_tax,
                               name = paste0(RANK_LABEL[[level]], " (within clade)")) +
    ggplot2::scale_colour_manual(values = cols, name = "Clade") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = if (by == "sample") "Sample" else "Group",
                  y = "Relative abundance",
                  title = "nosZ clade-resolved composition",
                  subtitle = paste0("Clade segments outlined by clade colour; ",
                                    RANK_LABEL[[level]], " nested within each clade")) +
    theme_nosz(base_size = base_size, x_angle = if (by == "sample") 45 else 0)
  p
}

#' Clade abundance ratio (e.g. nosZ-II / nosZ-I)
#'
#' nosZ Clade II / Clade I 
#' " nosZ  N2O " log2  + 
#'  =  = 
#'
#' @param x A `noszref` object.
#' @param numerator,denominator Branch names for the ratio.
#' @param palette Palette for groups.
#' @return A ggplot object.
#' @export
plot_clade_ratio <- function(x, numerator = "nosZII", denominator = "nosZI",
                             palette = "nosz") {
  stopifnot(inherits(x, "noszref"))
  need <- c(numerator, denominator)
  miss <- setdiff(need, x$branches)
  if (length(miss) > 0) {
    nosz_abort("", paste(miss, collapse = ", "),
               "", paste(x$branches, collapse = ", "))
  }
  a <- merge(x$rpkm[[numerator]], x$rpkm[[denominator]],
             by = "SampleID", suffixes = c("_num", "_den"))
  if (nrow(a) == 0) nosz_abort(" RPKM ")
  a <- merge(a, x$design[, c("SampleID", "Group")], by = "SampleID", all.x = TRUE)
  a$ratio <- log2((a$RPKM_num + 1e-6) / (a$RPKM_den + 1e-6))
  a$SampleID <- factor(a$SampleID, levels = order_samples(a$SampleID, x$design))
  cols <- group_palette(a$Group, palette)
  gmean <- stats::aggregate(ratio ~ Group, data = a, FUN = mean)

  ggplot2::ggplot(a, ggplot2::aes(SampleID, ratio, fill = Group)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
    ggplot2::geom_col(width = 0.65, colour = "grey25", linewidth = 0.2) +
    ggplot2::geom_hline(data = gmean, ggplot2::aes(yintercept = ratio, colour = Group),
                        linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(x = NULL,
                  y = expression(log[2] * " ratio"),
                  title = paste0(numerator, " / ", denominator, " abundance ratio"),
                  subtitle = paste0("dashed = group mean; positive values: ",
                                    numerator, " dominant")) +
    theme_nosz(x_angle = 45)
}

#' Per-sample clade statistics table
#'
#'  RPKM
#'
#' @param x A `noszref` object.
#' @return A data frame with `SampleID`, `Group`, per-branch RPKM and
#'   relative proportion columns.
#' @export
calc_clade_stats <- function(x) {
  stopifnot(inherits(x, "noszref"))
  d <- x$design[, c("SampleID", "Group"), drop = FALSE]
  for (b in x$branches) {
    r <- x$rpkm[[b]]
    if (is.null(r)) next
    d[[paste0(b, "_RPKM")]] <- r$RPKM[match(d$SampleID, r$SampleID)]
  }
  rk <- paste0(x$branches, "_RPKM")
  tot <- rowSums(d[, rk, drop = FALSE], na.rm = TRUE)
  for (b in x$branches) {
    v <- d[[paste0(b, "_RPKM")]]
    d[[paste0(b, "_prop")]] <- ifelse(tot > 0, v / tot, NA_real_)
  }
  d$Total_RPKM <- tot
  d
}

#' Pairwise clade abundance scatter
#'
#'  RPKM log10 Spearman 
#' 
#'
#' @param x A `noszref` object.
#' @param palette Palette for groups.
#' @return A ggplot (patchwork when >2 branches) object.
#' @export
plot_clade_scatter <- function(x, palette = "nosz") {
  stopifnot(inherits(x, "noszref"))
  bs <- x$branches[x$branches %in% names(x$rpkm)]
  if (length(bs) < 2) nosz_abort(" RPKM ")
  pairs <- utils::combn(bs, 2, simplify = FALSE)
  cols <- group_palette(x$design$Group, palette)
  ps <- lapply(pairs, function(pr) {
    a <- merge(x$rpkm[[pr[1]]], x$rpkm[[pr[2]]], by = "SampleID",
               suffixes = c(".x", ".y"))
    a <- merge(a, x$design[, c("SampleID", "Group")], by = "SampleID", all.x = TRUE)
    rho <- stats::cor.test(a$RPKM.x, a$RPKM.y, method = "spearman",
                           exact = FALSE)
    lab <- sprintf("rho = %.2f, P = %s", rho$estimate, fmt_p(rho$p.value))
    ggplot2::ggplot(a, ggplot2::aes(RPKM.x + 1e-6, RPKM.y + 1e-6, colour = Group)) +
      ggplot2::geom_point(size = 1.8, alpha = 0.9) +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(values = cols) +
      ggplot2::labs(x = paste0(pr[1], " RPKM (log10)"),
                    y = paste0(pr[2], " RPKM (log10)"),
                    subtitle = lab) +
      theme_nosz(legend_position = "none")
  })
  if (length(ps) == 1) return(ps[[1]])
  patchwork::wrap_plots(ps, ncol = length(ps)) &
    ggplot2::theme(legend.position = "bottom")
}
