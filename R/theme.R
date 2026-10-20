# ---------------------------------------------------------------------------
# 统一主题、调色板与导出
# ---------------------------------------------------------------------------

#' Built-in colour palettes of nosZviz
#'
#' @param n Number of colours.
#' @param palette Palette name, one of `"nosz"`, `"nosz_soft"`, `"Paired"`,
#'   `"Set3"`, `"npg"`, `"lancet"`, `"jama"`, `"tableau"`, `"viridis"`,
#'   `"earth"`, `"microbiome"`.
#' @return A character vector of colours.
#' @export
nosz_pal <- function(n = 3, palette = "nosz") {
  manual <- list(
    nosz = c("#1F78B4", "#FF7F00", "#33A02C", "#E31A1C", "#6A3D9A",
             "#A6CEE3", "#B15928", "#FB9A99", "#CAB2D6", "#FFFF99"),
    nosz_soft = c("#7AA8C4", "#E7A85C", "#8CC08A", "#D98C8C", "#9B87B8",
                  "#C9BE9E", "#8BBFC0", "#C7A0A0", "#A8B8C8", "#D8CAA0"),
    tableau = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
                "#86BCB6", "#D37295", "#A0CBE8", "#FFBE7D", "#8CD17D",
                "#F1CE63", "#D4A6C8", "#FABFD2", "#B6992D", "#499894"),
    earth = c("#3B5A6B", "#6C8A5E", "#A8B545", "#D9A441", "#C46A4A",
              "#8C5A7A", "#5F8A8B", "#B08A5A", "#7A8C6C", "#C7B8A0"),
    microbiome = c("#00A087", "#3C5488", "#E64B35", "#4DBBD5", "#F39B7F",
                   "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
  )
  brewer <- c("Paired", "Set3", "Set1", "Set2", "Dark2", "Accent", "Pastel1",
              "Pastel2", "Spectral", "RdYlBu")
  if (palette %in% names(manual)) {
    base <- manual[[palette]]
  } else if (palette %in% brewer) {
    k <- min(max(n, 3), 12)
    base <- RColorBrewer::brewer.pal(k, palette)
  } else if (identical(palette, "viridis")) {
    return(viridisLite::viridis(n))
  } else if (identical(palette, "npg")) {
    return(rep(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                 "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"),
               length.out = n))
  } else if (identical(palette, "lancet")) {
    return(rep(c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F",
                 "#FDAF91", "#AD002A", "#ADB6B6", "#1B1919"),
               length.out = n))
  } else if (identical(palette, "jama")) {
    return(rep(c("#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97",
                 "#6A6599", "#80796B"), length.out = n))
  } else {
    base <- manual$nosz
  }
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

#' Theme of nosZviz
#'
#' A compact, journal-ready theme: white background, thin black panel border,
#' no grid lines, outward ticks.
#'
#' @param base_size Base font size (mm-based figures usually use 7--9).
#' @param base_family Font family, default `"sans"`.
#' @param border Logical, draw panel border.
#' @param grid Logical, draw light grid lines.
#' @param x_angle Rotation angle of x axis labels.
#' @param legend_position Legend position.
#' @return A ggplot2 theme object.
#' @export
theme_nosz <- function(base_size = 8, base_family = "sans", border = TRUE,
                       grid = FALSE, x_angle = 0,
                       legend_position = "right") {
  th <- ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = if (grid) ggplot2::element_line(colour = "grey92", linewidth = 0.25) else ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = if (border) ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5) else ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size),
      axis.text.x = ggplot2::element_text(angle = x_angle,
                                          hjust = if (x_angle == 0) 0.5 else 1,
                                          vjust = if (x_angle == 0) 0.5 else 1),
      axis.title = ggplot2::element_text(colour = "black", size = base_size + 0.5),
      plot.title = ggplot2::element_text(size = base_size + 1.2, face = "bold",
                                         hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, colour = "grey25"),
      plot.tag = ggplot2::element_text(size = base_size + 1.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "black",
                                               linewidth = 0.4),
      strip.text = ggplot2::element_text(size = base_size, colour = "black"),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(size = base_size, face = "plain"),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key.size = grid::unit(3.2, "mm"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank()
    )
  th
}

#' Discrete fill scale with nosZviz palettes
#'
#' @param palette Palette name, see [nosz_pal()].
#' @param ... Passed to [ggplot2::discrete_scale()].
#' @return A ggplot2 scale.
#' @export
scale_fill_nosz <- function(palette = "nosz", ...) {
  ggplot2::discrete_scale("fill", pal = function(n) nosz_pal(n, palette), ...)
}

#' Discrete colour scale with nosZviz palettes
#'
#' @inheritParams scale_fill_nosz
#' @return A ggplot2 scale.
#' @export
scale_colour_nosz <- function(palette = "nosz", ...) {
  ggplot2::discrete_scale("colour", pal = function(n) nosz_pal(n, palette), ...)
}

# ggplot2 4.x 把 colour 拆成 colour/color 两个 aesthetic，这里补一个别名
#' @rdname scale_colour_nosz
#' @export
scale_color_nosz <- scale_colour_nosz

#' Save a figure in publication format
#'
#' Writes a figure to one or several formats with sizes expressed in
#' millimetres (default single-column 89 mm, double-column 183 mm).
#'
#' @param plot A ggplot / patchwork / grob object.
#' @param file Output file name without extension.
#' @param width,height Size in millimetres.
#' @param formats Character vector of formats: any of `"pdf"`, `"png"`,
#'   `"tiff"`, `"svg"`.
#' @param dpi Raster resolution.
#' @param bg Background colour.
#' @return Invisibly, the vector of written file paths.
#' @export
nosz_save <- function(plot, file, width = 89, height = 59,
                      formats = c("pdf", "png"), dpi = 600, bg = "white") {
  dir <- dirname(file)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  written <- character(0)
  for (fmt in formats) {
    path <- paste0(file, ".", fmt)
    # 绘制并保证设备一定被关闭（即使 print 出错），避免污染后续绘图
    ok <- tryCatch({
      if (fmt == "pdf") {
        grDevices::pdf(path, width = width / 25.4, height = height / 25.4,
                       useDingbats = FALSE, bg = bg)
        print(plot)
      } else if (fmt %in% c("png", "tiff", "jpeg")) {
        if (requireNamespace("ragg", quietly = TRUE) && fmt == "png") {
          ragg::agg_png(path, width = width / 25.4, height = height / 25.4,
                        units = "in", res = dpi, background = bg)
        } else {
          grDevices::png(path, width = width / 25.4, height = height / 25.4,
                         units = "in", res = dpi, bg = bg, type = "cairo")
        }
        print(plot)
      } else if (fmt == "svg") {
        grDevices::svg(path, width = width / 25.4, height = height / 25.4)
        print(plot)
      } else {
        next
      }
      TRUE
    }, error = function(e) {
      attr(e, "nosz_msg") <- conditionMessage(e)
      e
    }, finally = {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    })
    if (inherits(ok, "error")) {
      nosz_info("  [warn] ", basename(path), " : ",
                gsub("\n", " ", conditionMessage(ok)))
      next
    }
    if (isTRUE(ok)) written <- c(written, path)
  }
  if (length(written) == 0 && length(formats) > 0) {
    nosz_abort("图保存失败：", basename(file))
  }
  invisible(written)
}

# 配色便捷函数：按样本所属分组着色
sample_colours <- function(design, palette = "nosz", group_col = "Group") {
  if (is.null(design)) return(NULL)
  g <- design[[group_col]]
  pal <- group_palette(g, palette)
  stats::setNames(pal[as.character(g)], as.character(design$SampleID))
}
