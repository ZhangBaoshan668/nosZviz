# ---------------------------------------------------------------------------
#  NosZRef 
# ---------------------------------------------------------------------------

# beta otutab 
read_table_rownames <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  if (file.info(file)$size == 0) return(NULL)
  out <- try(
    utils::read.table(file, header = TRUE, row.names = 1, sep = "\t",
                      comment.char = "", quote = "", check.names = FALSE,
                      stringsAsFactors = FALSE, fill = TRUE),
    silent = TRUE
  )
  if (inherits(out, "try-error")) return(NULL)
  out <- as.matrix(out)
  storage.mode(out) <- "double"
  out[is.na(out)] <- 0
  out
}

#  SampleID/Group 
read_nosz_design <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  df <- read_delim_safe(file, delim = "\t")
  if (is.null(df) || ncol(df) < 2) {
    df <- read_delim_safe(file, delim = ",")
  }
  if (is.null(df) || ncol(df) < 2) return(NULL)
  nm <- tolower(colnames(df))
  id_i <- which(nm %in% c("sampleid", "sample", "sample_id", "id"))[1]
  gr_i <- which(nm %in% c("group", "groups", "treatment", "type"))[1]
  if (is.na(id_i)) id_i <- 1
  if (is.na(gr_i)) gr_i <- setdiff(seq_len(ncol(df)), id_i)[1]
  out <- data.frame(
    SampleID = as.character(df[[id_i]]),
    Group = as.character(df[[gr_i]]),
    stringsAsFactors = FALSE
  )
  extra <- setdiff(seq_len(ncol(df)), c(id_i, gr_i))
  for (i in extra) out[[colnames(df)[i]]] <- df[[i]]
  out <- out[!is.na(out$SampleID) & nzchar(out$SampleID), , drop = FALSE]
  out$Group[is.na(out$Group) | !nzchar(out$Group)] <- "All"
  out
}

#  merge  merge 
find_merge_dir <- function(path) {
  cand <- c(file.path(path, "merge"), path)
  hit <- cand[dir.exists(cand) & length(list.dirs(cand, recursive = FALSE)) > 0]
  hit <- hit[vapply(hit, function(p) {
    any(grepl("^nosZ", basename(list.dirs(p, recursive = FALSE))))
  }, logical(1))]
  if (length(hit) == 0) {
    nosz_abort(" `", path, "`  NosZRef  merge  nosZI/nosZII/nosZIII")
  }
  hit[1]
}

# 
discover_branches <- function(merge_dir) {
  dirs <- basename(list.dirs(merge_dir, recursive = FALSE))
  found <- dirs[grepl("^nosZ(I|II|III)$", dirs)]
  rank <- match(found, c("nosZI", "nosZII", "nosZIII"))
  found[order(rank)]
}

#  merge 
read_pe_sample_tables <- function(outdir, branches) {
  samp_dirs <- list.dirs(outdir, recursive = FALSE)
  samp_dirs <- samp_dirs[basename(samp_dirs) != "merge"]
  if (length(samp_dirs) == 0) return(NULL)
  out <- list()
  for (b in branches) {
    rp <- lapply(samp_dirs, function(d) {
      f <- file.path(d, paste0(b, "_", basename(d), ".rpkm.txt"))
      if (!file.exists(f)) return(NULL)
      x <- read_delim_safe(f, delim = "\t")
      if (is.null(x)) return(NULL)
      data.frame(SampleID = basename(d), RPKM = as.numeric(x[[2]]),
                 stringsAsFactors = FALSE)
    })
    rp <- do.call(rbind, rp)
    if (!is.null(rp) && nrow(rp) > 0) out[[b]] <- rp
  }
  out
}

# 
NOSZ_RANKS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
# 
RANK_LABEL <- c(kingdom = "Kingdom", phylum = "Phylum", class = "Class",
                order = "Order", family = "Family", genus = "Genus",
                species = "Species")

# rpkm / count
read_tax_level <- function(bdir, b, level, value = c("rpkm", "count")) {
  value <- match.arg(value)
  if (value == "rpkm") {
    cand <- c(file.path(bdir, paste0(b, "_merged_", level, "_rpkm.csv")),
              file.path(bdir, paste0(b, "_", level, "_merged_rpkm.csv")))
  } else {
    cand <- c(file.path(bdir, paste0(b, "_", level, "_merged_count.csv")),
              file.path(bdir, paste0(b, "_merged_", level, "_count.csv")))
  }
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) return(NULL)
  df <- read_delim_safe(f, delim = ",")
  if (is.null(df)) return(NULL)
  m <- df_to_matrix(df, 1)
  if (is.null(m) || ncol(m) == 0) return(NULL)
  m
}

# 
read_network <- function(ndir, b) {
  if (!dir.exists(ndir)) return(NULL)
  nodes_f <- file.path(ndir, "All_nodes_topology_attributes.csv")
  edges_f <- file.path(ndir, "All_edges_attributes.csv")
  stats_f <- file.path(ndir, "Network_summary_statistics.csv")
  if (!file.exists(nodes_f)) {
    per <- list.files(ndir, pattern = "^Node_attributes_.*csv$", full.names = TRUE)
    if (length(per) == 0) return(NULL)
    nodes_f <- per
  }
  nodes <- do.call(rbind, lapply(nodes_f, function(f) read_delim_safe(f, delim = ",")))
  edges <- if (file.exists(edges_f)) read_delim_safe(edges_f, delim = ",") else {
    per <- list.files(ndir, pattern = "^Edge_attributes_.*csv$", full.names = TRUE)
    if (length(per) > 0) {
      do.call(rbind, lapply(per, function(f) read_delim_safe(f, delim = ",")))
    } else NULL
  }
  stats <- if (file.exists(stats_f)) read_delim_safe(stats_f, delim = ",") else NULL
  if (is.null(nodes) || nrow(nodes) == 0) return(NULL)
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  colnames(nodes)[1] <- "Group"
  graph <- NULL
  if (!is.null(edges) && nrow(edges) > 0) {
    edges <- as.data.frame(edges, stringsAsFactors = FALSE)
    graph <- try(igraph::graph_from_data_frame(
      edges[, c("Source", "Target")], directed = FALSE,
      vertices = unique(nodes[, c("Group", "Node")])
    ), silent = TRUE)
    if (inherits(graph, "try-error")) graph <- NULL
  }
  list(nodes = nodes, edges = edges, stats = stats, graph = graph)
}

#' Read a NosZRef output directory
#'
#' @description
#' `read_noszref()` walks the whole result folder produced by the NosZRef
#' pipeline and collects every table that is needed for downstream
#' visualization: branch-level RPKM, read counts, OTU tables, rarefied OTU
#' tables, seven-rank taxonomic profiles (RPKM and counts), alpha diversity
#' indices, rarefaction curves, six beta diversity distance matrices, BLAST
#' taxonomic lineages and co-occurrence network topology.
#'
#' Missing steps are tolerated: a branch or an analysis that was not run is
#' simply absent, and plotting functions raise an informative error only when
#' they are actually asked for it.
#'
#' @param path Path to the NosZRef output directory (the one that contains
#'   `merge/`), or directly to the `merge/` directory.
#' @param design Sample metadata. Either a data frame with columns
#'   `SampleID` and `Group`, or a path to a two-column tab-delimited metadata
#'   file. If `NULL`, `<path>/metadata.txt` is used when present, otherwise all
#'   samples are assigned to one group called `"All"`.
#' @param branches Character vector of branches to read, e.g.
#'   `c("nosZI", "nosZII")`. Defaults to every branch found on disk.
#' @param network Logical, read co-occurrence network topology as well.
#' @param verbose Logical, print a short progress message.
#' @return An object of class `noszref`: a list with elements `path`,
#'   `design`, `branches`, `total_reads`, `rpkm`, `counts`, `otutab`,
#'   `otutab_rare`, `alpha`, `rare`, `beta`, `tax`, `tax_long`, `network`.
#' @export
#' @examples
#' \dontrun{
#' res <- read_noszref("path/to/NosZRef/Example_data")
#' res
#' }
read_noszref <- function(path, design = NULL, branches = NULL,
                         network = TRUE, verbose = TRUE) {
  if (!dir.exists(path)) nosz_abort("", path)
  merge_dir <- find_merge_dir(path)
  # merge  NosZRef metadata.txt 
  #  if(...) outdir <- merge_dir else outdir <- merge_dir 
  #    merge_dir metadata.txt 
  outdir <- if (identical(basename(merge_dir), "merge")) dirname(merge_dir) else merge_dir
  avail <- discover_branches(merge_dir)
  if (length(avail) == 0) nosz_abort(" `", merge_dir, "`  nosZ ")
  if (is.null(branches)) branches <- avail
  branches <- intersect(avail, branches)
  if (length(branches) == 0) {
    nosz_abort("", paste(avail, collapse = ", "))
  }

  # ----  ----
  # design_src  print() 
  # "" All
  design_src <- NULL
  if (is.null(design)) {
    auto <- c(file.path(outdir, "metadata.txt"), file.path(outdir, "group.txt"),
              file.path(outdir, "design.txt"), file.path(path, "metadata.txt"))
    auto <- unique(auto[file.exists(auto)])
    if (length(auto) > 0) {
      design_src <- auto[1]
      design <- read_nosz_design(design_src)
    } else {
      design <- NULL
    }
  } else if (is.character(design)) {
    design_src <- design
    design <- read_nosz_design(design)
  } else if (is.data.frame(design)) {
    design_src <- "<data.frame>"
    nm <- tolower(colnames(design))
    id_i <- which(nm %in% c("sampleid", "sample", "sample_id", "id"))[1]
    if (is.na(id_i)) id_i <- 1
    if (!"Group" %in% colnames(design)) {
      gr_i <- which(nm %in% c("group", "groups", "treatment"))[1]
      design$Group <- if (!is.na(gr_i)) design[[gr_i]] else "All"
    }
    design$SampleID <- as.character(design[[id_i]])
    design <- design[, c("SampleID", "Group", setdiff(colnames(design), c("SampleID", "Group")))]
  }

  # ----  reads ----
  tr_f <- file.path(merge_dir, "sample_merged_total_reads.txt")
  total_reads <- NULL
  if (file.exists(tr_f)) {
    x <- read_delim_safe(tr_f, delim = "\t", col_names = c("SampleID", "TotalReads"))
    if (!is.null(x)) total_reads <- x
  }
  if (is.null(total_reads)) {
    fs <- list.files(outdir, pattern = "_total_reads\\.txt$", recursive = TRUE,
                     full.names = TRUE)
    if (length(fs) > 0) {
      total_reads <- do.call(rbind, lapply(fs, function(f) {
        x <- read_delim_safe(f, delim = "\t", col_names = c("SampleID", "TotalReads"))
        if (is.null(x)) return(NULL)
        data.frame(SampleID = sub("_total_reads\\.txt$", "", basename(f)),
                   TotalReads = as.numeric(x[[2]]), stringsAsFactors = FALSE)
      }))
    }
  }
  if (!is.null(total_reads)) {
    total_reads <- as.data.frame(total_reads, stringsAsFactors = FALSE)
    total_reads$SampleID <- as.character(total_reads$SampleID)
    total_reads$TotalReads <- as.numeric(total_reads$TotalReads)
  }

  # ----  design ----
  if (is.null(design) || nrow(design) == 0) {
    smp <- NULL
    if (!is.null(total_reads)) smp <- total_reads$SampleID
    if (is.null(smp)) {
      f <- file.path(merge_dir, branches[1], paste0(branches[1], "_merged_rpkm.txt"))
      if (file.exists(f)) smp <- read_delim_safe(f, delim = "\t")[[1]]
    }
    design <- data.frame(SampleID = as.character(smp), Group = "All",
                         stringsAsFactors = FALSE)
    # 
    # ""
    if (is.null(design_src)) {
      warning(" ",
              paste0(c(file.path(outdir, "metadata.txt"), file.path(outdir, "group.txt"),
                       file.path(outdir, "design.txt"), file.path(path, "metadata.txt")),
                     collapse = ""),
              " ", nrow(design), "  \"All\"",
              "LEfSebeta ",
              " read_noszref(..., design = \"metadata.txt\") ",
              call. = FALSE)
    }
  }

  # ----  ----
  res <- list(
    path = normalizePath(merge_dir, mustWork = FALSE),
    outdir = normalizePath(outdir, mustWork = FALSE),
    design = design,
    design_file = design_src,
    branches = branches,
    total_reads = total_reads,
    rpkm = list(), counts = list(), otutab = list(), otutab_rare = list(),
    alpha = list(), rare = list(), beta = list(), tax = list(),
    tax_long = list(), network = list()
  )

  for (b in branches) {
    bdir <- file.path(merge_dir, b)

    #  RPKM /  reads
    f <- file.path(bdir, paste0(b, "_merged_rpkm.txt"))
    if (file.exists(f)) {
      x <- read_delim_safe(f, delim = "\t")
      if (!is.null(x)) {
        res$rpkm[[b]] <- data.frame(SampleID = as.character(x[[1]]),
                                    RPKM = as.numeric(x[[2]]),
                                    stringsAsFactors = FALSE)
      }
    }
    f <- file.path(bdir, paste0(b, "_merged_results.txt"))
    if (file.exists(f)) {
      x <- read_delim_safe(f, delim = "\t", col_names = c("SampleID", "Count"))
      if (!is.null(x)) {
        res$counts[[b]] <- data.frame(SampleID = as.character(x[[1]]),
                                      Count = as.numeric(x[[2]]),
                                      stringsAsFactors = FALSE)
      }
    }

    # OTU 
    f <- file.path(bdir, paste0(b, "_otutab.txt"))
    m <- read_table_rownames(f)
    if (!is.null(m)) res$otutab[[b]] <- m

    #  alpha
    f <- file.path(bdir, "alpha", paste0(b, "_otutab_rare.txt"))
    m <- read_table_rownames(f)
    if (!is.null(m)) res$otutab_rare[[b]] <- m

    f <- file.path(bdir, "alpha", paste0(b, "_alpha_div.txt"))
    if (file.exists(f)) {
      x <- read_delim_safe(f, delim = "\t")
      if (!is.null(x)) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        colnames(x)[1] <- "SampleID"
        x$SampleID <- as.character(x$SampleID)
        for (j in setdiff(colnames(x), "SampleID")) x[[j]] <- as.numeric(x[[j]])
        res$alpha[[b]] <- x
      }
    }

    f <- file.path(bdir, "alpha", paste0(b, "_alpha_rare.txt"))
    if (file.exists(f)) {
      x <- read_delim_safe(f, delim = "\t")
      if (!is.null(x)) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        colnames(x)[1] <- "Depth"
        x$Depth <- as.numeric(x$Depth)
        res$rare[[b]] <- x
      }
    }

    # beta 
    bdir_beta <- file.path(bdir, "beta")
    if (dir.exists(bdir_beta)) {
      files <- list.files(bdir_beta, pattern = "\\.txt$", full.names = TRUE)
      files <- files[!grepl("sorted", files)]
      mats <- list()
      for (f in files) {
        nm <- sub("\\.txt$", "", basename(f))
        m <- read_table_rownames(f)
        if (!is.null(m) && nrow(m) > 1) mats[[nm]] <- m
      }
      if (length(mats) > 0) res$beta[[b]] <- mats
    }

    # 
    lv <- list()
    for (level in NOSZ_RANKS) {
      r <- read_tax_level(bdir, b, level, "rpkm")
      c0 <- read_tax_level(bdir, b, level, "count")
      if (!is.null(r) || !is.null(c0)) lv[[level]] <- list(rpkm = r, count = c0)
    }
    if (length(lv) > 0) res$tax[[b]] <- lv

    # BLAST 
    f <- file.path(bdir, paste0(b, "_merged_tax.txt"))
    if (file.exists(f) && file.info(f)$size > 0) {
      x <- try(readr::read_delim(
        f, delim = "\t", col_names = FALSE, show_col_types = FALSE,
        progress = FALSE, quote = "", trim_ws = TRUE,
        col_types = readr::cols(.default = readr::col_character())
      ), silent = TRUE)
      if (!inherits(x, "try-error") && ncol(x) >= 15) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        nm <- c("qseqid", "qlen", "lineage", "sgi", "slen", "pident", "alen",
                "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore", "staxid", "ssciname")
        colnames(x)[seq_len(min(length(nm), ncol(x)))] <- nm[seq_len(min(length(nm), ncol(x)))]
        x$pident <- as.numeric(x$pident)
        x$bitscore <- as.numeric(x$bitscore)
        #  qseqid  design 
        ids <- unique(design$SampleID)
        ids <- ids[order(nchar(ids), decreasing = TRUE)]
        x$SampleID <- vapply(x$qseqid, function(q) {
          hit <- ids[vapply(ids, function(i) startsWith(q, paste0(i, "_")), logical(1))]
          if (length(hit) == 0) NA_character_ else hit[1]
        }, character(1))
        res$tax_long[[b]] <- x
      }
    }

    # 
    if (network) {
      nw <- read_network(file.path(bdir, "network"), b)
      if (!is.null(nw)) res$network[[b]] <- nw
    }
  }

  #  merge  rpkm
  if (length(res$rpkm) == 0) {
    pe <- read_pe_sample_tables(outdir, branches)
    if (!is.null(pe)) res$rpkm <- pe
  }

  #  ID design 
  # /
  data_ids <- unique(unlist(lapply(res$rpkm, `[[`, "SampleID")))
  if (length(data_ids) > 0) {
    miss <- setdiff(res$design$SampleID, data_ids)
    if (length(miss) > 0) {
      warning("design  ", length(miss), " ",
              paste(utils::head(miss, 8), collapse = ", "),
              if (length(miss) > 8) " ..." else "",
              " rpkm ", call. = FALSE)
    }
    extra <- setdiff(data_ids, res$design$SampleID)
    if (length(extra) > 0) {
      warning(" ", length(extra), "  design ",
              paste(utils::head(extra, 8), collapse = ", "),
              if (length(extra) > 8) " ..." else "", call. = FALSE)
    }
    # design  metadata 
  }

  structure(res, class = "noszref")
}

#' @export
print.noszref <- function(x, ...) {
  cat("NosZRef result object\n")
  cat("  path     :", x$path, "\n")
  cat("  branches :", paste(x$branches, collapse = ", "), "\n")
  grp <- unique(x$design$Group)
  cat("  samples  : ", nrow(x$design), " (", paste(grp, collapse = "/"), ")\n", sep = "")
  cat("  design   :", if (is.null(x$design_file)) "<none, all samples in \"All\">" else x$design_file, "\n")
  for (b in x$branches) {
    cat("  --", b, "\n")
    cat("      rpkm:", ifelse(is.null(x$rpkm[[b]]), "no", "yes"),
        " otutab:", ifelse(is.null(x$otutab[[b]]), "no", "yes"),
        " alpha:", ifelse(is.null(x$alpha[[b]]), "no", "yes"),
        " beta:", ifelse(is.null(x$beta[[b]]), "no", length(x$beta[[b]])),
        " rank:", ifelse(is.null(x$tax[[b]]), "no", length(x$tax[[b]])),
        " network:", ifelse(is.null(x$network[[b]]), "no", "yes"), "\n")
  }
  invisible(x)
}

#' Overview table of a noszref object
#'
#' @param x A `noszref` object.
#' @return A data frame summarizing availability of each result type.
#' @export
noszref_summary <- function(x) {
  stopifnot(inherits(x, "noszref"))
  out <- lapply(x$branches, function(b) {
    data.frame(
      Branch = b,
      Samples = if (!is.null(x$rpkm[[b]])) nrow(x$rpkm[[b]]) else NA_integer_,
      TotalRPKM = if (!is.null(x$rpkm[[b]])) sum(x$rpkm[[b]]$RPKM, na.rm = TRUE) else NA_real_,
      OTU = if (!is.null(x$otutab[[b]])) nrow(x$otutab[[b]]) else NA_integer_,
      TaxonomyRanks = if (!is.null(x$tax[[b]])) length(x$tax[[b]]) else 0L,
      BetaDistances = if (!is.null(x$beta[[b]])) length(x$beta[[b]]) else 0L,
      NetworkNodes = if (!is.null(x$network[[b]])) nrow(x$network[[b]]$nodes) else 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Extract a taxon x sample matrix from a noszref object
#'
#' @param x A `noszref` object.
#' @param branch Branch name, e.g. `"nosZI"`.
#' @param level Taxonomic level, one of `"kingdom"`, `"phylum"`, `"class"`,
#'   `"order"`, `"family"`, `"genus"`, `"species"`.
#' @param value `"rpkm"` (default) or `"count"`.
#' @return A numeric matrix with taxa in rows and samples in columns.
#' @export
nosz_tax_matrix <- function(x, branch, level = "genus", value = "rpkm") {
  stopifnot(inherits(x, "noszref"))
  tx <- nosz_need(x, "tax", branch)
  if (is.null(tx[[level]])) {
    nosz_abort(" ", branch, "  ", level, " ",
               paste(names(tx), collapse = ", "))
  }
  m <- tx[[level]][[value]]
  if (is.null(m)) m <- tx[[level]][["rpkm"]]
  if (is.null(m)) nosz_abort(" ", branch, "  ", level, " ")
  #  design 
  smp <- intersect(x$design$SampleID, colnames(m))
  m[, c(smp, setdiff(colnames(m), smp)), drop = FALSE]
}

#' OTU table of one branch
#'
#' @inheritParams nosz_tax_matrix
#' @param rarefied Logical, return the rarefied table when available.
#' @return A numeric matrix (OTU x sample).
#' @export
nosz_otu_matrix <- function(x, branch, rarefied = FALSE) {
  stopifnot(inherits(x, "noszref"))
  m <- if (rarefied) x$otutab_rare[[branch]] else x$otutab[[branch]]
  if (is.null(m)) {
    nosz_abort(" ", branch, " ", ifelse(rarefied, "", ""), " OTU ")
  }
  smp <- intersect(x$design$SampleID, colnames(m))
  m[, c(smp, setdiff(colnames(m), smp)), drop = FALSE]
}
