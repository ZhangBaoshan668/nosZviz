# ---------------------------------------------------------------------------
# 模块六：共现网络
# ---------------------------------------------------------------------------

# 由节点/边表重建 igraph 对象
network_graph <- function(x, branch, group = NULL) {
  nw <- nosz_need(x, "network", branch)
  nodes <- nw$nodes
  edges <- nw$edges
  if (is.null(edges) || nrow(edges) == 0) nosz_abort("分支 ", branch, " 没有边表。")
  if (!is.null(group)) {
    nodes <- nodes[nodes$Group == group, , drop = FALSE]
    edges <- edges[edges$Group == group, , drop = FALSE]
  }
  verts <- data.frame(name = as.character(nodes$Node), nodes,
                      stringsAsFactors = FALSE)
  verts <- verts[!duplicated(verts$name), , drop = FALSE]
  e <- data.frame(from = as.character(edges$Source),
                  to = as.character(edges$Target),
                  Cor = suppressWarnings(as.numeric(edges$Correlation)),
                  W = suppressWarnings(as.numeric(edges$Weight)),
                  stringsAsFactors = FALSE)
  e <- e[e$from %in% verts$name & e$to %in% verts$name, , drop = FALSE]
  e <- e[e$from != e$to, , drop = FALSE]
  g <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = verts)
  g <- igraph::simplify(g, edge.attr.comb = list(Cor = "mean", W = "mean", "ignore"))
  g
}

# 节点属性整理：模块取前 N 个，其余归为 Other
node_attrs <- function(g, top_modules = 8) {
  v <- igraph::as_data_frame(g, what = "vertices")
  v$Degree <- suppressWarnings(as.numeric(v$Degree))
  v$Degree[is.na(v$Degree)] <- igraph::degree(g)
  v$Modularity <- as.character(v$Modularity)
  tb <- sort(table(v$Modularity), decreasing = TRUE)
  keep <- names(tb)[seq_len(min(top_modules, length(tb)))]
  v$Module <- ifelse(v$Modularity %in% keep,
                     paste0("M", v$Modularity), "Other")
  v$Module <- factor(v$Module, levels = c(paste0("M", keep), "Other"))
  v$Hub <- !is.na(v$Role) & v$Role != "Peripheral"
  v
}

#' Co-occurrence network plot
#'
#' Draws the group-wise co-occurrence network reconstructed from the
#' `Node_attributes_*.csv` and `Edge_attributes_*.csv` files written by
#' NosZRef. Nodes are coloured by module, sized by degree, and hub nodes
#' (Zi > 2.5 or Pi > 0.62) are highlighted and labelled.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param group One group, or `NULL` for the first group in the design.
#' @param layout igraph layout name, e.g. `"fr"`, `"kk"`, `"stress"`,
#'   `"circle"`.
#' @param colour_by `"module"` or `"role"`.
#' @param node_size Node size scaling factor.
#' @param label_hubs Logical, label hub nodes.
#' @param label_all Logical, label every node (only for small networks).
#' @param top_modules Number of modules kept as separate colours.
#' @param seed Random seed for the layout.
#' @return A ggplot object.
#' @export
plot_network <- function(x, branch, group = NULL, layout = "fr",
                         colour_by = c("module", "role"), node_size = 1,
                         label_hubs = TRUE, label_all = FALSE,
                         top_modules = 8, seed = 42) {
  stopifnot(inherits(x, "noszref"))
  colour_by <- match.arg(colour_by)
  nw <- nosz_need(x, "network", branch)
  groups <- unique(nw$nodes$Group)
  group <- group %||% groups[1]
  if (!group %in% groups) nosz_abort("分组 `", group, "` 不在网络结果中。可用：",
                                     paste(groups, collapse = ", "))
  g <- network_graph(x, branch, group)
  v <- node_attrs(g, top_modules)
  igraph::V(g)$Module <- v$Module[match(igraph::V(g)$name, v$name)]
  igraph::V(g)$Hub <- v$Hub[match(igraph::V(g)$name, v$name)]
  igraph::V(g)$Role <- v$Role[match(igraph::V(g)$name, v$name)]
  igraph::V(g)$Degree2 <- v$Degree[match(igraph::V(g)$name, v$name)]
  ed <- igraph::as_data_frame(g, what = "edges")
  igraph::E(g)$Sign <- ifelse(!is.na(ed$Cor) & ed$Cor > 0, "Positive",
                              ifelse(!is.na(ed$Cor) & ed$Cor < 0, "Negative", "NA"))

  set.seed(seed)
  lay <- ggraph::create_layout(g, layout = layout)
  # 只给实际出现的模块配色，避免 levels 与取值数量不一致
  mods <- unique(as.character(v$Module))
  n_real <- sum(mods != "Other")
  fill_cols <- stats::setNames(rep("grey82", length(mods)), mods)
  if (n_real > 0) {
    fill_cols[mods != "Other"] <- nosz_pal(n_real, "tableau")[seq_len(n_real)]
  }
  role_cols <- c(`Network hub` = "#D7191C", `Module hub` = "#FDAE61",
                 Connector = "#2C7AB6", Peripheral = "grey80")

  p <- ggraph::ggraph(lay) +
    ggraph::geom_edge_link(ggplot2::aes(edge_colour = Sign),
                           linewidth = 0.12, alpha = 0.28) +
    ggraph::scale_edge_colour_manual(
      values = c(Positive = "#D73027", Negative = "#4575B4", `NA` = "grey75"),
      name = "Correlation", guide = ggplot2::guide_legend(override.aes = list(alpha = 1))
    )
  if (colour_by == "module") {
    p <- p + ggraph::geom_node_point(ggplot2::aes(size = Degree2, fill = Module),
                                     shape = 21, colour = "grey25", stroke = 0.18) +
      ggplot2::scale_fill_manual(values = fill_cols, name = "Module")
  } else {
    p <- p + ggraph::geom_node_point(ggplot2::aes(size = Degree2, fill = Role),
                                     shape = 21, colour = "grey25", stroke = 0.18) +
      ggplot2::scale_fill_manual(values = role_cols, name = "Role")
  }
  if (label_hubs) {
    p <- p + ggraph::geom_node_text(ggplot2::aes(filter = Hub, label = name),
                                    size = 1.55, repel = TRUE, max.overlaps = 25,
                                    seed = 1, colour = "black")
  }
  if (label_all) {
    p <- p + ggraph::geom_node_text(ggplot2::aes(label = name), size = 1.1,
                                    repel = TRUE, max.overlaps = 50, seed = 1)
  }
  st <- nw$stats
  sub <- if (!is.null(st) && group %in% st$Group) {
    r <- st[st$Group == group, , drop = FALSE]
    paste0("Nodes = ", r$Nodes, ", Edges = ", r$Edges,
           ", Modularity = ", formatC(r$Modularity, 3, format = "f"))
  } else {
    paste0("Nodes = ", igraph::vcount(g), ", Edges = ", igraph::ecount(g))
  }
  p +
    ggplot2::scale_size_continuous(range = c(0.4, 2.6) * node_size,
                                   name = "Degree") +
    ggplot2::labs(x = NULL, y = NULL, subtitle = sub,
                  title = paste0(branch, " co-occurrence network - ", group)) +
    theme_nosz(border = FALSE, legend_position = "right") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   legend.key.size = grid::unit(2.6, "mm"))
}

#' Zi-Pi plot for node topological roles
#'
#' The standard module-hub classification: nodes with Zi > 2.5 are module hubs,
#' Pi > 0.62 are connectors, both are network hubs.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param label_n Number of nodes to label (highest Zi + Pi).
#' @param zi_cut,pi_cut Thresholds.
#' @param palette Palette name.
#' @return A ggplot object.
#' @export
plot_network_zipi <- function(x, branch, label_n = 10, zi_cut = 2.5,
                              pi_cut = 0.62, palette = "nosz") {
  stopifnot(inherits(x, "noszref"))
  nw <- nosz_need(x, "network", branch)
  v <- nw$nodes
  v$Zi <- suppressWarnings(as.numeric(v$Zi))
  v$Pi <- suppressWarnings(as.numeric(v$Pi))
  v <- v[!is.na(v$Zi) & !is.na(v$Pi), , drop = FALSE]
  if (nrow(v) == 0) nosz_abort("网络节点缺少 Zi/Pi 值。")
  v$Role <- factor(v$Role, levels = c("Network hub", "Module hub", "Connector", "Peripheral"))
  role_cols <- c(`Network hub` = "#D7191C", `Module hub` = "#FDAE61",
                 Connector = "#2C7AB6", Peripheral = "grey75")
  lab <- v[order(-(v$Zi + v$Pi)), , drop = FALSE]
  lab <- head(lab[lab$Role != "Peripheral", , drop = FALSE], label_n)
  ggplot2::ggplot(v, ggplot2::aes(Zi, Pi)) +
    ggplot2::geom_hline(yintercept = pi_cut, linetype = "dashed",
                        linewidth = 0.25, colour = "grey45") +
    ggplot2::geom_vline(xintercept = zi_cut, linetype = "dashed",
                        linewidth = 0.25, colour = "grey45") +
    ggplot2::geom_point(ggplot2::aes(colour = Role, size = Degree), alpha = 0.85) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = Node), size = 1.6,
                             max.overlaps = 25, seed = 1,
                             min.segment.length = 0.1, segment.size = 0.15) +
    ggplot2::scale_colour_manual(values = role_cols, name = NULL,
                                 drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(0.3, 2.2), name = "Degree") +
    ggplot2::annotate("text", x = zi_cut, y = max(v$Pi, pi_cut) * 1.02,
                      label = "Module hub", hjust = -0.05, size = 1.9,
                      colour = "grey30") +
    ggplot2::annotate("text", x = min(v$Zi), y = pi_cut,
                      label = "Connector", vjust = -0.5, hjust = 0, size = 1.9,
                      colour = "grey30") +
    ggplot2::facet_wrap(~Group) +
    ggplot2::labs(x = "Zi (within-module connectivity)",
                  y = "Pi (among-module connectivity)",
                  title = paste0(branch, " - node topological roles")) +
    theme_nosz(grid = TRUE, legend_position = "right")
}

#' Network topological properties
#'
#' Compares nodes, edges, average degree, modularity, transitivity, diameter
#' and average path length among groups.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param metrics Metrics to show; `NULL` for all.
#' @param palette Palette name.
#' @param ncol Panel columns.
#' @return A ggplot object.
#' @export
plot_network_properties <- function(x, branch, metrics = NULL,
                                    palette = "nosz", ncol = 4) {
  stopifnot(inherits(x, "noszref"))
  nw <- nosz_need(x, "network", branch)
  st <- nw$stats
  if (is.null(st)) nosz_abort("缺少 Network_summary_statistics.csv。")
  metrics <- metrics %||% intersect(c("Nodes", "Edges", "Average_degree",
                                      "Modularity", "Transitivity", "Diameter",
                                      "Average_path_length"), colnames(st))
  long <- tidyr::pivot_longer(st, dplyr::all_of(metrics), names_to = "Metric",
                              values_to = "Value")
  long$Metric <- factor(long$Metric, levels = metrics,
                        labels = gsub("_", " ", metrics))
  long$Group <- factor(long$Group, levels = unique(st$Group))
  ggplot2::ggplot(long, ggplot2::aes(Group, Value, fill = Group)) +
    ggplot2::geom_col(width = 0.65, colour = "black", linewidth = 0.22) +
    ggplot2::geom_text(ggplot2::aes(label = formatC(Value, digits = 2,
                                                    format = "f")),
                       vjust = -0.35, size = 1.7) +
    ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = group_palette(long$Group, palette)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0(branch, " - network topological properties")) +
    theme_nosz(x_angle = 30, legend_position = "none")
}

#' Degree distribution with power-law fit
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param palette Palette name.
#' @param fit Logical, add a log-log linear fit per group.
#' @return A ggplot object.
#' @export
plot_network_degree <- function(x, branch, palette = "nosz", fit = TRUE) {
  stopifnot(inherits(x, "noszref"))
  nw <- nosz_need(x, "network", branch)
  v <- nw$nodes
  v$Degree <- suppressWarnings(as.numeric(v$Degree))
  v <- v[!is.na(v$Degree) & v$Degree > 0, , drop = FALSE]
  if (nrow(v) == 0) nosz_abort("没有可用的度值。")
  tbl <- v %>% dplyr::group_by(Group, Degree) %>%
    dplyr::summarise(N = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(P = N / sum(N))
  cols <- group_palette(tbl$Group, palette)
  p <- ggplot2::ggplot(tbl, ggplot2::aes(Degree, P, colour = Group)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::annotation_logticks(sides = "bl", linewidth = 0.2) +
    ggplot2::labs(x = "Degree (k)", y = "P(k)",
                  title = paste0(branch, " - degree distribution")) +
    theme_nosz(grid = TRUE, legend_position = "top")
  if (fit) {
    fits <- tbl %>% dplyr::filter(Degree > 1, P > 0) %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(
        m = stats::coef(stats::lm(log10(P) ~ log10(Degree)))[2],
        r2 = summary(stats::lm(log10(P) ~ log10(Degree)))$r.squared,
        .groups = "drop"
      )
    p <- p + ggplot2::geom_smooth(ggplot2::aes(group = Group), method = "lm",
                                  formula = y ~ x, se = FALSE, linewidth = 0.35,
                                  alpha = 0.6)
    p <- p + ggplot2::labs(subtitle = paste0(
      paste0(fits$Group, ": \u03b3 = ", formatC(fits$m, 2, format = "f"),
             " (R\u00b2 = ", formatC(fits$r2, 2, format = "f"), ")"),
      collapse = "; "))
  }
  p
}

#' Module size and connectivity
#'
#' Compares the size and mean degree of the largest modules between groups.
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param top_modules Number of modules shown per group.
#' @param palette Palette name.
#' @return A ggplot object.
#' @export
plot_network_modules <- function(x, branch, top_modules = 8,
                                 palette = "tableau") {
  stopifnot(inherits(x, "noszref"))
  nw <- nosz_need(x, "network", branch)
  v <- nw$nodes
  v$Degree <- suppressWarnings(as.numeric(v$Degree))
  v$Modularity <- as.character(v$Modularity)
  agg <- v %>% dplyr::group_by(Group, Modularity) %>%
    dplyr::summarise(Nodes = dplyr::n(),
                     MeanDegree = mean(Degree, na.rm = TRUE), .groups = "drop")
  agg <- agg %>% dplyr::group_by(Group) %>%
    dplyr::arrange(dplyr::desc(Nodes), .by_group = TRUE) %>%
    dplyr::mutate(rank = dplyr::row_number(),
                  Module = ifelse(rank <= top_modules,
                                  paste0("M", Modularity), "Other")) %>%
    dplyr::ungroup()
  agg$Module <- factor(agg$Module,
                       levels = c(paste0("M", sort(unique(v$Modularity))),
                                  "Other"))
  long <- tidyr::pivot_longer(agg, c("Nodes", "MeanDegree"),
                              names_to = "Metric", values_to = "Value")
  long$Metric <- factor(long$Metric, labels = c("Module size (nodes)",
                                                "Mean degree within module"))
  n_mod <- length(levels(droplevels(agg$Module)))
  cols <- stats::setNames(c(nosz_pal(max(n_mod - 1, 1), palette), "grey82"),
                          levels(droplevels(agg$Module)))
  ggplot2::ggplot(long, ggplot2::aes(Module, Value, fill = Module)) +
    ggplot2::geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
    ggplot2::facet_grid(Group ~ Metric, scales = "free_y") +
    ggplot2::scale_fill_manual(values = cols, name = "Module") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0(branch, " - module size and connectivity")) +
    theme_nosz(x_angle = 45, legend_position = "none")
}

#' UpSet plot of shared network nodes between groups
#'
#' @param x A `noszref` object.
#' @param branch Branch name.
#' @param top_n Number of intersections shown.
#' @param palette Bar colour.
#' @return A patchwork object.
#' @export
plot_network_upset <- function(x, branch, top_n = 15,
                               palette = "#4E79A7") {
  stopifnot(inherits(x, "noszref"))
  nw <- nosz_need(x, "network", branch)
  nodes <- nw$nodes
  sets <- split(as.character(nodes$Node), as.character(nodes$Group))
  sets <- sets[vapply(sets, length, integer(1)) > 0]
  if (length(sets) < 2) nosz_abort("至少需要两个分组的网络才能做交集分析。")
  alln <- unique(unlist(sets, use.names = FALSE))
  memb <- vapply(sets, function(s) alln %in% s, logical(length(alln)))
  if (is.null(dim(memb))) memb <- matrix(memb, ncol = length(sets))
  colnames(memb) <- names(sets)
  key <- apply(memb, 1, function(r) paste(colnames(memb)[r], collapse = " & "))
  tb <- as.data.frame(table(key), stringsAsFactors = FALSE)
  colnames(tb) <- c("Combo", "Freq")
  tb <- tb[order(-tb$Freq), , drop = FALSE]
  if (nrow(tb) > top_n) {
    rest <- sum(tb$Freq[(top_n + 1):nrow(tb)])
    tb <- rbind(tb[seq_len(top_n), , drop = FALSE],
                data.frame(Combo = "Others", Freq = rest))
  }
  tb$Combo <- factor(tb$Combo, levels = tb$Combo)
  # 点阵数据
  grid <- expand.grid(Combo = levels(tb$Combo), Set = names(sets),
                      stringsAsFactors = FALSE)
  grid$In <- vapply(seq_len(nrow(grid)), function(i) {
    cb <- as.character(grid$Combo[i]); st <- grid$Set[i]
    if (cb == "Others") return(NA)
    st %in% strsplit(cb, " & ", fixed = TRUE)[[1]]
  }, logical(1))
  grid$Combo <- factor(grid$Combo, levels = levels(tb$Combo))
  grid$Set <- factor(grid$Set, levels = rev(names(sets)))

  p_top <- ggplot2::ggplot(tb, ggplot2::aes(Combo, Freq)) +
    ggplot2::geom_col(width = 0.65, fill = palette, colour = "black",
                      linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = Freq), vjust = -0.35, size = 1.7) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = "Shared nodes") +
    theme_nosz(x_angle = 45, legend_position = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  p_bot <- ggplot2::ggplot(grid, ggplot2::aes(Combo, Set)) +
    ggplot2::geom_point(data = grid[is.na(grid$In), , drop = FALSE],
                        colour = "grey85", size = 0.8) +
    ggplot2::geom_point(data = grid[!is.na(grid$In) & !grid$In, , drop = FALSE],
                        colour = "grey85", size = 0.8) +
    ggplot2::geom_point(data = grid[!is.na(grid$In) & grid$In, , drop = FALSE],
                        ggplot2::aes(colour = Set), size = 1.5) +
    ggplot2::geom_line(data = grid[!is.na(grid$In) & grid$In, , drop = FALSE],
                       ggplot2::aes(group = Combo), colour = "grey35",
                       linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = group_palette(names(sets), "nosz"),
                                 guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nosz(x_angle = 45, legend_position = "none", grid = TRUE)

  patchwork::wrap_plots(p_top, p_bot, ncol = 1, heights = c(1.6, 1)) +
    patchwork::plot_annotation(
      title = paste0(branch, " - shared and unique network nodes"),
      subtitle = paste0(paste(names(sets), collapse = " / "),
                        "; total ", length(alln), " unique nodes")
    )
}
