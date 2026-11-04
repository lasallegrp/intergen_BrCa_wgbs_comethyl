#!/usr/bin/env Rscript

# ============================================================
# 12B_enrich_selected_modules.R
# ============================================================
#
# ENRICHMENT METHODS AVAILABLE
#   - KEGG        : clusterProfiler::enrichKEGG()    — offline-ish
#   - GO BP/MF/CC : clusterProfiler::enrichGO()      — fully offline
#   - Reactome    : ReactomePA::enrichPathway()       — fully offline
#   - Enrichr     : enrichR::enrichr()               — requires internet
#
# SELECTION MODES
#   Provide either:
#     --module_list_file  : explicit list of module names to enrich
#                           (can be a single .txt file OR a directory of .txt files;
#                            if a directory, each .txt is processed separately and
#                            outputs are placed in a subfolder named after the file stem)
#   OR:
#     --trait_list_file + --stats_file_v1 (+ v2/v3)
#       modules selected where trait in trait_list AND p <= p_thresh
#       (can also be a single .txt file OR a directory of .txt files;
#        same output-subfolder behaviour applies)
#
# REQUIRED INPUTS
#   --project_root      : root directory of the project
#   --annotation_dir_v1 : path to v1 annotation folder from Script 12a
#
# OPTIONAL INPUTS
#   --annotation_dir_v2       : v2 annotation folder
#   --annotation_dir_v3       : v3 annotation folder
#   --module_list_file        : text file (or directory of .txt files) with module names
#   --stats_file_v1/v2/v3     : ME-trait stats TSV from Script 09a
#   --trait_list_file         : text file (or directory of .txt files) with trait names
#   --p_thresh                : module selection p threshold [default = 0.05]
#   --do_kegg                 : TRUE/FALSE [default = true]
#   --do_go                   : TRUE/FALSE [default = true]
#   --do_reactome             : TRUE/FALSE [default = true]
#   --do_enrichr              : TRUE/FALSE [default = false]
#   --enrichr_dbs             : comma-separated Enrichr DB names
#   --organism                : KEGG organism code [default = hsa]
#   --use_enrichr_background  : TRUE/FALSE [default = true]
#   --min_bg_genes            : minimum background size [default = 200]
#   --enrichr_retries         : Enrichr retry count [default = 3]
#   --enrichr_sleep           : seconds between retries [default = 1]
#   --plot_p_col              : p.adjust or pvalue for plots [default = p.adjust]
#   --plot_top_n              : top N terms per module in summary plots [default = 10]
#   --plot_sig_cutoff         : significance cutoff for summary plots [default = 0.05]
#   --plot_width              : summary plot width in inches [default = 14]
#   --plot_height             : summary plot height in inches [default = 10]
#
# OUTPUTS
#   project_root/comethyl_output/12b_enrichment/<cpg>/<region>/<variant>/[<list_stem>/]
#     selected_modules.txt
#     selected_modules_summary.tsv
#     <module>_KEGG.tsv / .xlsx / _dotplot.pdf
#     <module>_GO_BP/MF/CC.tsv / .xlsx / _dotplot.pdf
#     <module>_Reactome.tsv / .xlsx / _dotplot.pdf
#     <module>_enrichr.xlsx
#     <module>_enrichr_<DB>.pdf          <- per-module per-database plotEnrich PDFs
#     Summary_KEGG_top<n>_<pcol><cutoff>.pdf
#     Summary_GO_BP/MF/CC_top<n>_<pcol><cutoff>.pdf
#     Summary_Reactome_top<n>_<pcol><cutoff>.pdf
#     Summary_Enrichr_<db>_top<n>_<pcol><cutoff>.pdf
#     run_log.txt
#     run_parameters.txt
# ============================================================

message("Starting Script 12b")

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(scales)
})

# ============================================================
# 1) Basic helpers
# ============================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx  <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

trim_or_null <- function(x) {
  if (is.null(x) || is.na(x)) return(NULL)
  x <- trimws(x)
  if (!nzchar(x)) return(NULL)
  x
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
}

stop_if_missing <- function(x, label) {
  if (is.null(x) || !nzchar(x)) stop("Missing required argument: ", label, call. = FALSE)
}

validate_file_exists <- function(path, label) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(as.character(path))) {
    stop(label, " is missing or invalid. Received: ",
         paste(capture.output(str(path)), collapse = " "), call. = FALSE)
  }
  path <- as.character(path)
  if (!file.exists(path) || dir.exists(path)) {
    stop(label, " not found as a file: ", path, call. = FALSE)
  }
  invisible(normalizePath(path, mustWork = TRUE))
}

validate_dir_exists <- function(path, label) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(as.character(path))) {
    stop(label, " is missing or invalid. Received: ",
         paste(capture.output(str(path)), collapse = " "), call. = FALSE)
  }
  path <- as.character(path)
  if (!dir.exists(path)) {
    stop(label, " not found as a directory: ", path, call. = FALSE)
  }
  invisible(normalizePath(path, mustWork = TRUE))
}

timestamp_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

append_log <- function(logfile, ...) {
  txt <- paste0("[", timestamp_now(), "] ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = logfile, append = TRUE)
}

write_lines_safe <- function(x, file) {
  writeLines(as.character(x), con = file, useBytes = TRUE)
}

qc_gene_list <- function(x, label = "genes") {
  x0 <- x
  x  <- trimws(as.character(x))
  x  <- x[!is.na(x) & x != ""]
  x  <- unique(x)
  ok <- grepl("^[A-Za-z0-9._-]+$", x)
  message(sprintf("[%s] input n=%d -> cleaned unique n=%d", label, length(x0), length(x)))
  message(sprintf("[%s] invalid-format count=%d", label, sum(!ok)))
  list(clean = x, ok = ok)
}

# ============================================================
# 2a) File-or-directory resolver
#
#   resolve_list_input(path) returns a named character vector:
#     names  = file stem (used as output subfolder name, or "" for single-file mode)
#     values = absolute path to the .txt file
#
#   Rules:
#     - NULL input              -> returns NULL (caller handles "not provided")
#     - path is a single .txt   -> returns c("" = path)  (no extra subfolder)
#     - path is a directory     -> returns one entry per *.txt in that dir,
#                                  named by tools::file_path_sans_ext(basename(...))
# ============================================================
resolve_list_input <- function(path) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(as.character(path))) {
    return(NULL)
  }

  path <- as.character(path)

  if (file.exists(path) && !dir.exists(path)) {
    # Single-file mode: no additional output subfolder is created.
    return(setNames(normalizePath(path, mustWork = TRUE), ""))
  }

  if (dir.exists(path)) {
    # Directory mode: each .txt file is processed separately.
    # Output is written to a subfolder named after the text-file stem.
    txts <- sort(list.files(path, pattern = "\\.txt$", full.names = TRUE))
    if (length(txts) == 0)
      stop("Directory provided for list file contains no .txt files: ", path, call. = FALSE)
    stems <- tools::file_path_sans_ext(basename(txts))
    return(setNames(normalizePath(txts, mustWork = TRUE), stems))
  }

  stop("Path does not exist as a file or directory: ", path, call. = FALSE)
}

# ============================================================
# 2b) Input readers
# ============================================================
read_clean_txt_list <- function(path) {
  validate_file_exists(path, basename(path))
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!grepl("^#", x)]   # allow comment lines in list files
  unique(x)
}

read_module_list_file <- function(path) {
  read_clean_txt_list(path)
}

read_trait_list_file <- function(path) {
  read_clean_txt_list(path)
}

load_annotation_files <- function(annotation_dir) {
  annotated_regions_file <- file.path(annotation_dir, "Annotated_Regions.tsv")
  gene_list_file         <- file.path(annotation_dir, "Module_Gene_List.tsv")
  module_summary_file    <- file.path(annotation_dir, "Module_Gene_Summary.tsv")
  background_file        <- file.path(annotation_dir, "Background_Genes.tsv")

  validate_file_exists(annotated_regions_file, "Annotated_Regions.tsv")
  validate_file_exists(gene_list_file,         "Module_Gene_List.tsv")
  validate_file_exists(module_summary_file,    "Module_Gene_Summary.tsv")

  annotated_regions <- read.delim(annotated_regions_file,
                                  stringsAsFactors = FALSE, check.names = FALSE)
  gene_list         <- read.delim(gene_list_file,
                                  stringsAsFactors = FALSE, check.names = FALSE)
  module_summary    <- read.delim(module_summary_file,
                                  stringsAsFactors = FALSE, check.names = FALSE)

  background_genes <- NULL
  if (file.exists(background_file)) {
    bg <- read.delim(background_file, stringsAsFactors = FALSE, check.names = FALSE)
    if ("gene_symbol" %in% names(bg)) background_genes <- bg$gene_symbol
  }

  miss_regions <- setdiff(c("RegionID", "chr", "start", "end", "module"),
                          colnames(annotated_regions))
  if (length(miss_regions) > 0)
    stop("Annotated_Regions.tsv missing required columns: ",
         paste(miss_regions, collapse = ", "), call. = FALSE)

  miss_gene <- setdiff(c("module", "gene_symbol"), colnames(gene_list))
  if (length(miss_gene) > 0)
    stop("Module_Gene_List.tsv missing required columns: ",
         paste(miss_gene, collapse = ", "), call. = FALSE)

  list(annotated_regions = annotated_regions,
       gene_list         = gene_list,
       module_summary    = module_summary,
       background_genes  = background_genes)
}

# ============================================================
# 3) Output directory helper
# ============================================================
derive_pipeline_dirs_from_annotation_dir <- function(annotation_dir,
                                                      project_root,
                                                      step_name,
                                                      list_stem = "") {
  variant_name  <- basename(annotation_dir)
  region_label  <- basename(dirname(annotation_dir))
  cpg_label     <- basename(dirname(dirname(annotation_dir)))

  pipeline_root <- file.path(project_root, "comethyl_output")
  step_dir      <- file.path(pipeline_root, step_name)

  if (nzchar(list_stem)) {
    out_dir <- file.path(step_dir, cpg_label, region_label, variant_name, list_stem)
  } else {
    out_dir <- file.path(step_dir, cpg_label, region_label, variant_name)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  list(pipeline_root = pipeline_root, step_dir = step_dir, out_dir = out_dir,
       variant_name = variant_name, region_label = region_label,
       cpg_label = cpg_label)
}

# ============================================================
# 4) Stats-based module selection
# ============================================================
standardize_stats_columns <- function(df) {
  cn <- tolower(colnames(df))

  for (col in c("module", "trait")) {
    if (!col %in% cn) stop("stats file must contain ", col, " column.", call. = FALSE)
    names(df)[match(col, cn)] <- col
    cn <- tolower(colnames(df))
  }

  p_candidates <- c("p", "pvalue", "p.value", "p_value")
  found_p      <- p_candidates[p_candidates %in% cn][1]
  if (is.na(found_p))
    stop("stats file must contain a p-value column named: ",
         paste(p_candidates, collapse = ", "), call. = FALSE)
  names(df)[match(found_p, tolower(colnames(df)))] <- "p"
  df
}

read_stats_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df  <- switch(ext,
    xlsx = openxlsx::read.xlsx(path),
    xls  = openxlsx::read.xlsx(path),
    tsv  = ,
    txt  = read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    csv  = read.csv(path,   stringsAsFactors = FALSE, check.names = FALSE),
    stop("Unsupported stats file format: ", path, call. = FALSE)
  )
  standardize_stats_columns(df)
}

select_modules_from_stats <- function(stats_df, trait_list, p_thresh) {
  stats_df %>%
    dplyr::mutate(p = as.numeric(p)) %>%
    dplyr::filter(!is.na(p), trait %in% trait_list, p <= p_thresh) %>%
    dplyr::pull(module) %>%
    unique() %>%
    as.character()
}

# ============================================================
# 5) Gene list builder
# ============================================================
build_gene_vector_for_module <- function(gene_list_df, module_name) {
  gene_list_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene_symbol) %>%
    as.character() %>%
    trimws() %>%
    unique() %>%
    .[!is.na(.) & . != ""]
}

# ============================================================
# 6) KEGG enrichment — offline-ish via clusterProfiler
# ============================================================
run_kegg_clusterprofiler <- function(gene_symbols, out_prefix, organism = "hsa") {
  gene_symbols <- unique(na.omit(as.character(gene_symbols)))
  gene_symbols <- gene_symbols[nzchar(gene_symbols)]
  if (length(gene_symbols) == 0) return(invisible(NULL))

  gene_df <- suppressMessages(
    clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(gene_df) || nrow(gene_df) == 0) return(invisible(NULL))

  kegg_res <- suppressMessages(
    clusterProfiler::enrichKEGG(gene = unique(gene_df$ENTREZID),
                                 organism = organism,
                                 pvalueCutoff = 1, qvalueCutoff = 1)
  )
  kegg_df <- as.data.frame(kegg_res)
  if (is.null(kegg_df) || nrow(kegg_df) == 0) return(invisible(NULL))

  write.table(kegg_df, paste0(out_prefix, "_KEGG.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  openxlsx::write.xlsx(kegg_df, paste0(out_prefix, "_KEGG.xlsx"), rowNames = FALSE)

  top <- kegg_df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(25, nrow(kegg_df)))
  if (nrow(top) > 0) {
    top$Description <- factor(top$Description, levels = rev(top$Description))
    p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
      geom_point() + theme_bw(base_size = 12) +
      labs(title = "Top KEGG pathways (clusterProfiler)",
           x = "-log10(adj. p)", y = NULL)
    ggsave(paste0(out_prefix, "_KEGG_dotplot.pdf"), plot = p, width = 9, height = 6)
  }

  message("[KEGG] Done: ", nrow(kegg_df), " pathways saved")
  invisible(kegg_df)
}

# ============================================================
# 7) GO enrichment (BP, MF, CC) — fully offline via clusterProfiler
# ============================================================
run_go_clusterprofiler <- function(gene_symbols, out_prefix) {
  gene_symbols <- unique(na.omit(as.character(gene_symbols)))
  gene_symbols <- gene_symbols[nzchar(gene_symbols)]
  if (length(gene_symbols) == 0) return(invisible(NULL))

  gene_df <- suppressMessages(
    clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    message("[GO] No Entrez IDs mapped for: ", basename(out_prefix))
    return(invisible(NULL))
  }

  entrez_ids <- unique(gene_df$ENTREZID)

  for (ont in c("BP", "MF", "CC")) {
    res <- tryCatch({
      suppressMessages(
        clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db,
                                   ont = ont, pvalueCutoff = 1, qvalueCutoff = 1,
                                   readable = TRUE)
      )
    }, error = function(e) {
      message("[GO ", ont, "] Failed: ", conditionMessage(e)); NULL
    })

    if (is.null(res)) next
    df <- as.data.frame(res)
    if (nrow(df) == 0) next

    write.table(df, paste0(out_prefix, "_GO_", ont, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    openxlsx::write.xlsx(df, paste0(out_prefix, "_GO_", ont, ".xlsx"), rowNames = FALSE)

    top <- df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(25, nrow(df)))
    if (nrow(top) > 0) {
      top$Description <- factor(top$Description, levels = rev(top$Description))
      p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
        geom_point() + theme_bw(base_size = 12) +
        labs(title = paste0("Top GO ", ont, " terms (clusterProfiler)"),
             x = "-log10(adj. p)", y = NULL)
      ggsave(paste0(out_prefix, "_GO_", ont, "_dotplot.pdf"), plot = p, width = 9, height = 6)
    }
    message("[GO ", ont, "] Done: ", nrow(df), " terms saved")
  }
  invisible(NULL)
}

# ============================================================
# 8) Reactome enrichment — fully offline via ReactomePA
# ============================================================
run_reactome_clusterprofiler <- function(gene_symbols, out_prefix) {
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    message("[Reactome] ReactomePA not installed; skipping. ",
            "Install with: BiocManager::install('ReactomePA')")
    return(invisible(NULL))
  }

  gene_symbols <- unique(na.omit(as.character(gene_symbols)))
  gene_symbols <- gene_symbols[nzchar(gene_symbols)]
  if (length(gene_symbols) == 0) return(invisible(NULL))

  gene_df <- suppressMessages(
    clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    message("[Reactome] No Entrez IDs mapped for: ", basename(out_prefix))
    return(invisible(NULL))
  }

  res <- tryCatch({
    suppressMessages(
      ReactomePA::enrichPathway(gene = unique(gene_df$ENTREZID),
                                organism = "human", pvalueCutoff = 1,
                                qvalueCutoff = 1, readable = TRUE)
    )
  }, error = function(e) {
    message("[Reactome] Failed: ", conditionMessage(e)); NULL
  })

  if (is.null(res)) return(invisible(NULL))
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(invisible(NULL))

  write.table(df, paste0(out_prefix, "_Reactome.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  openxlsx::write.xlsx(df, paste0(out_prefix, "_Reactome.xlsx"), rowNames = FALSE)

  top <- df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(25, nrow(df)))
  if (nrow(top) > 0) {
    top$Description <- factor(top$Description, levels = rev(top$Description))
    p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
      geom_point() + theme_bw(base_size = 12) +
      labs(title = "Top Reactome Pathways (ReactomePA)",
           x = "-log10(adj. p)", y = NULL)
    ggsave(paste0(out_prefix, "_Reactome_dotplot.pdf"), plot = p, width = 9, height = 6)
  }

  message("[Reactome] Done: ", nrow(df), " pathways saved")
  invisible(df)
}

# ============================================================
# 9) Summary dotplot helpers — shared style across all methods
# ============================================================
load_method_results <- function(out_dir, modules, file_suffix,
                                p_col     = "p.adjust",
                                count_col = "Count",
                                desc_col  = "Description") {
  all_data <- data.frame()

  for (m in modules) {
    f <- file.path(out_dir, paste0(m, file_suffix))
    if (!file.exists(f)) next

    df <- tryCatch(openxlsx::read.xlsx(f), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next

    actual_p <- if (p_col %in% colnames(df)) {
      p_col
    } else if ("p.adjust" %in% colnames(df)) {
      message("[load] '", p_col, "' not in ", basename(f), " — using 'p.adjust'")
      "p.adjust"
    } else if ("pvalue" %in% colnames(df)) {
      message("[load] '", p_col, "' not in ", basename(f), " — using 'pvalue'")
      "pvalue"
    } else {
      message("[load] No p-value column found in ", basename(f), " — skipping")
      next
    }

    needed <- c(desc_col, actual_p, count_col)
    if (!all(needed %in% colnames(df))) {
      message("[load] Missing columns in ", basename(f),
              ". Need: ", paste(needed, collapse = ", "))
      next
    }

    df <- df %>%
      dplyr::select(Description = !!desc_col,
                    p_col_val   = !!actual_p,
                    gene_count  = !!count_col) %>%
      dplyr::mutate(module   = m,
                    p_source = actual_p)

    all_data <- dplyr::bind_rows(all_data, df)
  }
  all_data
}

make_enrichment_dotplot <- function(df,
                                    title,
                                    source_tag  = NULL,
                                    p_col_label = "adj. p",
                                    sig_cutoff  = 0.05,
                                    top_n       = 10,
                                    max_chars   = 55,
                                    plot_width  = 14,
                                    plot_height = 10,
                                    base_size   = 13,
                                    out_file    = NULL) {
  if (nrow(df) == 0) {
    message("[plot] No data to plot for: ", title)
    return(invisible(NULL))
  }

  df$Description <- ifelse(
    nchar(df$Description) > max_chars,
    paste0(substr(df$Description, 1, max_chars - 3), "..."),
    df$Description
  )

  df <- df %>%
    dplyr::filter(!is.na(p_col_val), p_col_val <= sig_cutoff) %>%
    dplyr::group_by(module) %>%
    dplyr::arrange(p_col_val, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  if (nrow(df) == 0) {
    message("[plot] No terms pass sig_cutoff=", sig_cutoff, " for: ", title)
    return(invisible(NULL))
  }

  df <- df %>%
    dplyr::arrange(module, p_col_val) %>%
    dplyr::mutate(Description = factor(Description,
                                       levels = rev(unique(Description))))

  source_label  <- if (!is.null(source_tag)) paste0(" | Source: ", source_tag) else ""
  subtitle_text <- paste0("Top ", top_n, " terms per module | ",
                          p_col_label, " \u2264 ", sig_cutoff, source_label)

  p <- ggplot(df, aes(x = module, y = Description,
                       size = gene_count, fill = p_col_val)) +
    geom_point(shape = 23, color = "grey30") +
    scale_size(range = c(3, 10), name = "Gene Count") +
    scale_fill_gradient(
      low    = "#B22222",
      high   = "#FFC0CB",
      name   = p_col_label,
      limits = c(0, sig_cutoff),
      oob    = scales::squish
    ) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x        = element_text(angle = 30, hjust = 1),
      axis.text.y        = element_text(size = base_size - 1),
      axis.title         = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "right"
    ) +
    labs(title = title, subtitle = subtitle_text)

  if (!is.null(out_file)) {
    ggsave(out_file, plot = p, width = plot_width, height = plot_height,
           dpi = 300, units = "in")
    message("[plot] Saved: ", basename(out_file))
  }
  invisible(p)
}

plot_kegg_summary <- function(out_dir, modules,
                               p_col = "p.adjust", sig_cutoff = 0.05,
                               top_n = 10, plot_width = 14, plot_height = 10,
                               base_size = 13) {
  p_col_label <- if (p_col == "p.adjust") "adj. p" else "raw p"
  p_col_safe  <- gsub("\\.", "", p_col)
  df          <- load_method_results(out_dir, modules, "_KEGG.xlsx", p_col = p_col)
  out_file    <- file.path(out_dir,
    sprintf("Summary_KEGG_top%d_%s%s.pdf", top_n, p_col_safe, sig_cutoff))
  make_enrichment_dotplot(df, title = "Top KEGG Pathways per Module",
                          source_tag = "clusterProfiler", p_col_label = p_col_label,
                          sig_cutoff = sig_cutoff, top_n = top_n,
                          plot_width = plot_width, plot_height = plot_height,
                          base_size = base_size, out_file = out_file)
}

plot_go_summary <- function(out_dir, modules, ontologies = c("BP", "MF", "CC"),
                             p_col = "p.adjust", sig_cutoff = 0.05,
                             top_n = 10, plot_width = 14, plot_height = 10,
                             base_size = 13) {
  p_col_label <- if (p_col == "p.adjust") "adj. p" else "raw p"
  p_col_safe  <- gsub("\\.", "", p_col)
  for (ont in ontologies) {
    df <- load_method_results(out_dir, modules,
                              paste0("_GO_", ont, ".xlsx"), p_col = p_col)
    if (nrow(df) == 0) next
    out_file <- file.path(out_dir,
      sprintf("Summary_GO_%s_top%d_%s%s.pdf", ont, top_n, p_col_safe, sig_cutoff))
    make_enrichment_dotplot(df, title = paste0("Top GO ", ont, " Terms per Module"),
                            source_tag = "clusterProfiler", p_col_label = p_col_label,
                            sig_cutoff = sig_cutoff, top_n = top_n,
                            plot_width = plot_width, plot_height = plot_height,
                            base_size = base_size, out_file = out_file)
  }
}

plot_reactome_summary <- function(out_dir, modules,
                                   p_col = "p.adjust", sig_cutoff = 0.05,
                                   top_n = 10, plot_width = 14, plot_height = 10,
                                   base_size = 13) {
  p_col_label <- if (p_col == "p.adjust") "adj. p" else "raw p"
  p_col_safe  <- gsub("\\.", "", p_col)
  df          <- load_method_results(out_dir, modules, "_Reactome.xlsx", p_col = p_col)
  out_file    <- file.path(out_dir,
    sprintf("Summary_Reactome_top%d_%s%s.pdf", top_n, p_col_safe, sig_cutoff))
  make_enrichment_dotplot(df, title = "Top Reactome Pathways per Module",
                          source_tag = "clusterProfiler / ReactomePA",
                          p_col_label = p_col_label,
                          sig_cutoff = sig_cutoff, top_n = top_n,
                          plot_width = plot_width, plot_height = plot_height,
                          base_size = base_size, out_file = out_file)
}

plot_enrichr_summary <- function(out_dir, modules, databases,
                                  p_col = "p.adjust", sig_cutoff = 0.05,
                                  top_n = 10, plot_width = 14, plot_height = 10,
                                  base_size = 13) {
  enrichr_p_col <- if (p_col == "p.adjust") "Adjusted.P.value" else "P.value"
  p_col_label   <- if (p_col == "p.adjust") "adj. p" else "raw p"
  p_col_safe    <- gsub("\\.", "", p_col)

  for (db in databases) {
    all_data <- data.frame()

    for (m in modules) {
      f <- file.path(out_dir, paste0(m, "_enrichr.xlsx"))
      if (!file.exists(f)) next

      sheet <- substr(gsub("[^A-Za-z0-9]", "_", db), 1, 31)
      df <- tryCatch(openxlsx::read.xlsx(f, sheet = sheet),
                     error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) next

      actual_p <- if (enrichr_p_col %in% colnames(df)) enrichr_p_col else
                  if ("Adjusted.P.value" %in% colnames(df)) "Adjusted.P.value" else "P.value"

      if (!all(c("Term", actual_p, "Overlap") %in% colnames(df))) next

      df$gene_count  <- suppressWarnings(
        as.integer(sub("/.*$", "", as.character(df$Overlap)))
      )
      df$Description <- gsub(" \\(GO:\\d+\\)$", "", df$Term)

      df <- df %>%
        dplyr::select(Description, p_col_val = !!actual_p, gene_count) %>%
        dplyr::mutate(module = m)

      all_data <- dplyr::bind_rows(all_data, df)
    }

    if (nrow(all_data) == 0) {
      message("[plot_enrichr] No data for DB: ", db)
      next
    }

    db_safe  <- gsub("[^A-Za-z0-9]+", "_", db)
    out_file <- file.path(out_dir,
      sprintf("Summary_Enrichr_%s_top%d_%s%s.pdf", db_safe, top_n, p_col_safe, sig_cutoff))

    make_enrichment_dotplot(all_data,
                            title      = paste0("Top Enrichr Terms: ", db),
                            source_tag = "Enrichr / maayanlab.cloud",
                            p_col_label = p_col_label,
                            sig_cutoff = sig_cutoff, top_n = top_n,
                            plot_width = plot_width, plot_height = plot_height,
                            base_size  = base_size, out_file = out_file)
  }
}

# ============================================================
# 10) Enrichr helpers
# ============================================================

init_enrichr <- function() {
  options(
    enrichR.sites.base.address = "https://maayanlab.cloud/",
    enrichR.base.address       = "https://maayanlab.cloud/Enrichr/",
    speedrichr.base.address    = "https://maayanlab.cloud/speedrichr/api/",
    enrichR.live               = TRUE,
    enrichR.quiet              = TRUE,
    modEnrichR.use             = TRUE,
    enrichR.sites              = c("Enrichr", "FlyEnrichr", "WormEnrichr",
                                   "YeastEnrichr", "FishEnrichr", "OxEnrichr")
  )
}

validate_enrichr_dbs_safe <- function(dbs) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    warning("enrichR package not installed; cannot validate Enrichr DBs.")
    return(list(ok = dbs, bad = character(0), avail = character(0)))
  }

  out <- tryCatch({
    init_enrichr()
    Sys.sleep(0.5)
    avail       <- enrichR::listEnrichrDbs()
    avail_names <- as.character(avail$libraryName)
    list(ok    = intersect(dbs, avail_names),
         bad   = setdiff(dbs, avail_names),
         avail = avail_names)
  }, error = function(e) {
    warning("Could not validate Enrichr DBs; using requested list as-is. Reason: ",
            conditionMessage(e))
    list(ok = dbs, bad = character(0), avail = character(0))
  })
  out
}

run_enrichr_simple <- function(gene_symbols, out_prefix, dbs,
                               background_genes = NULL,
                               retries          = 3,
                               sleep_time       = 1) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    warning("enrichR package not installed; skipping Enrichr.")
    return(invisible(NULL))
  }

  fg_qc        <- qc_gene_list(gene_symbols, "foreground")
  gene_symbols <- fg_qc$clean
  if (length(gene_symbols) == 0) return(invisible(NULL))
  if (length(dbs) == 0)          return(invisible(NULL))

  if (!is.null(background_genes)) {
    bg_qc            <- qc_gene_list(background_genes, "background")
    background_genes <- bg_qc$clean
    if (length(background_genes) == 0) background_genes <- NULL
  }

  init_enrichr()

  res <- NULL; last_err <- NULL

  for (attempt in seq_len(retries)) {
    message("[Enrichr] Attempt ", attempt, " of ", retries,
            " for: ", basename(out_prefix))

    res <- tryCatch({
      if (!is.null(background_genes) && length(background_genes) > 0) {
        enrichR::enrichr(genes = gene_symbols, databases = dbs,
                         background = background_genes)
      } else {
        enrichR::enrichr(genes = gene_symbols, databases = dbs)
      }
    }, error = function(e) { last_err <<- conditionMessage(e); NULL })

    if (!is.null(res) && length(res) > 0) break
    if (attempt < retries) {
      message("[Enrichr] Failed attempt ", attempt, ": ", last_err,
              " — retrying in ", sleep_time, "s")
      Sys.sleep(sleep_time)
    }
  }

  if (is.null(res) || length(res) == 0) {
    warning("Enrichr failed for ", out_prefix,
            if (!is.null(last_err)) paste0(": ", last_err) else ".")
    return(invisible(NULL))
  }

  # ---- Save combined Excel workbook ----
  wb <- openxlsx::createWorkbook(); wrote_any <- FALSE

  for (nm in names(res)) {
    df <- res[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    sheet <- substr(gsub("[^A-Za-z0-9]", "_", nm), 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df)
    wrote_any <- TRUE
  }

  if (wrote_any) {
    openxlsx::saveWorkbook(wb, paste0(out_prefix, "_enrichr.xlsx"), overwrite = TRUE)
    message("[Enrichr] Saved: ", basename(out_prefix), "_enrichr.xlsx")

    # ---- Per-database plotEnrich() PDFs ----
    # Requires enrichR::plotEnrich(); mirrors the old plot_and_save() logic.
    for (nm in names(res)) {
      df <- res[[nm]]
      if (is.null(df) || nrow(df) == 0) next

      db_safe  <- substr(gsub("[^A-Za-z0-9]", "_", nm), 1, 60)
      pdf_file <- paste0(out_prefix, "_enrichr_", db_safe, ".pdf")
      mod_name <- basename(out_prefix)   # e.g. "skyblue2"

      tryCatch({
        grDevices::pdf(pdf_file, height = 7, width = 15)
        print(
          enrichR::plotEnrich(df, showTerms = 25, numChar = 75,
                              y = "Count", orderBy = "P.value") +
            ggplot2::ggtitle(paste(nm, "for", mod_name, "module"))
        )
        grDevices::dev.off()
        message("[Enrichr] Plot saved: ", basename(pdf_file))
      }, error = function(e) {
        # Ensure graphics device is always closed on error
        tryCatch(grDevices::dev.off(), error = function(e2) NULL)
        message("[Enrichr] Plot failed for DB '", nm, "': ", conditionMessage(e))
      })
    }
  }

  Sys.sleep(sleep_time)
  invisible(res)
}

# ============================================================
# 11) Per-variant enrichment runner
# ============================================================
run_enrichment_for_variant <- function(annotation_dir,
                                       variant_label,
                                       module_list_file       = NULL,
                                       stats_file             = NULL,
                                       trait_list             = NULL,
                                       p_thresh               = 0.05,
                                       do_kegg                = TRUE,
                                       do_go                  = TRUE,
                                       do_reactome            = TRUE,
                                       do_enrichr             = FALSE,
                                       enrichr_dbs            = character(0),
                                       organism               = "hsa",
                                       project_root           = NULL,
                                       use_enrichr_background = TRUE,
                                       min_bg_genes           = 200,
                                       enrichr_retries        = 3,
                                       enrichr_sleep          = 1,
                                       plot_p_col             = "p.adjust",
                                       plot_top_n             = 10,
                                       plot_sig_cutoff        = 0.05,
                                       plot_width             = 14,
                                       plot_height            = 10,
                                       list_stem              = "") {

  validate_dir_exists(annotation_dir, paste0("annotation_dir_", variant_label))

  files            <- load_annotation_files(annotation_dir)
  gene_list_df     <- files$gene_list
  module_summary   <- files$module_summary
  background_genes <- files$background_genes

  dir_info <- derive_pipeline_dirs_from_annotation_dir(
    annotation_dir = annotation_dir,
    project_root   = project_root,
    step_name      = "12b_enrichment",
    list_stem      = list_stem
  )

  out_dir     <- dir_info$out_dir
  log_file    <- file.path(out_dir, "run_log.txt")
  params_file <- file.path(out_dir, "run_parameters.txt")
  selected_modules_file         <- file.path(out_dir, "selected_modules.txt")
  selected_modules_summary_file <- file.path(out_dir, "selected_modules_summary.tsv")

  append_log(log_file, "Starting enrichment for variant: ", variant_label,
             if (nzchar(list_stem)) paste0(" / list: ", list_stem) else "")
  append_log(log_file, "annotation_dir: ", annotation_dir)
  append_log(log_file, "out_dir: ", out_dir)
  append_log(log_file, "do_kegg: ", do_kegg, " | do_go: ", do_go,
             " | do_reactome: ", do_reactome, " | do_enrichr: ", do_enrichr)
  append_log(log_file, "plot_p_col: ", plot_p_col)

  # ---- Module selection ----
  selected_modules <- character(0)
  selection_mode   <- NULL

  if (!is.null(module_list_file)) {
    validate_file_exists(module_list_file, "module_list_file")
    selected_modules <- read_module_list_file(module_list_file)
    selection_mode   <- "module_list_file"
  } else if (!is.null(stats_file) && !is.null(trait_list) && length(trait_list) > 0) {
    validate_file_exists(stats_file, paste0("stats_file_", variant_label))
    stats_df         <- read_stats_file(stats_file)
    selected_modules <- select_modules_from_stats(stats_df, trait_list, p_thresh)
    selection_mode   <- "stats_file_plus_trait_list"
  } else {
    stop("For ", variant_label,
         ", provide --module_list_file OR both stats_file and trait_list_file.",
         call. = FALSE)
  }

  selected_modules <- unique(as.character(selected_modules))
  selected_modules <- selected_modules[nzchar(selected_modules)]
  selected_modules <- selected_modules[selected_modules %in% unique(gene_list_df$module)]
  selected_modules <- selected_modules[!grepl("^grey$", selected_modules, ignore.case = TRUE)]

  write_lines_safe(selected_modules, selected_modules_file)
  sel_summary <- module_summary %>% dplyr::filter(module %in% selected_modules)
  write.table(sel_summary, selected_modules_summary_file,
              sep = "\t", quote = FALSE, row.names = FALSE)

  append_log(log_file, "selection_mode: ", selection_mode)
  append_log(log_file, "selected modules count: ", length(selected_modules))
  append_log(log_file, "selected modules: ",
             ifelse(length(selected_modules) == 0, "(none)",
                    paste(selected_modules, collapse = ", ")))

  # ---- Enrichr background ----
  bg_to_use <- NULL; bg_ok <- FALSE

  if (isTRUE(use_enrichr_background) && !is.null(background_genes)) {
    bg_qc            <- qc_gene_list(background_genes, "background")
    background_genes <- bg_qc$clean

    mapped        <- AnnotationDbi::select(org.Hs.eg.db, keys = background_genes,
                                           keytype = "SYMBOL",
                                           columns = c("SYMBOL", "ENTREZID"))
    mapped_unique <- unique(mapped$SYMBOL[!is.na(mapped$ENTREZID)])
    map_rate      <- if (length(background_genes) > 0)
                       length(mapped_unique) / length(background_genes) else 0

    append_log(log_file, "background genes: ", length(background_genes))
    append_log(log_file, sprintf("background Entrez mapping rate: %.1f%%", 100 * map_rate))

    if (length(background_genes) >= min_bg_genes) {
      bg_to_use <- background_genes; bg_ok <- TRUE
    } else {
      append_log(log_file, "Background too small (", length(background_genes),
                 " < ", min_bg_genes, "); disabling background mode.")
    }
  } else {
    append_log(log_file, "Enrichr background mode disabled.")
  }

  if (length(selected_modules) == 0) {
    append_log(log_file, "No modules selected — skipping enrichment.")
    return(invisible(NULL))
  }

  # ---- Per-module enrichment ----
  for (m in selected_modules) {
    gene_vec   <- build_gene_vector_for_module(gene_list_df, m)
    out_prefix <- file.path(out_dir, m)

    append_log(log_file, "Module ", m, ": genes=", length(gene_vec))
    message("\n--- Enriching module: ", m, " (", length(gene_vec), " genes) ---")

    if (do_kegg)     { message("[", m, "] KEGG...");     run_kegg_clusterprofiler(gene_vec, out_prefix, organism) }
    if (do_go)       { message("[", m, "] GO...");       run_go_clusterprofiler(gene_vec, out_prefix) }
    if (do_reactome) { message("[", m, "] Reactome..."); run_reactome_clusterprofiler(gene_vec, out_prefix) }
    if (do_enrichr) {
      message("[", m, "] Enrichr...")
      run_enrichr_simple(gene_vec, out_prefix, enrichr_dbs,
                         background_genes = if (bg_ok) bg_to_use else NULL,
                         retries    = enrichr_retries,
                         sleep_time = enrichr_sleep)
    }
  }

  # ---- Summary plots across all selected modules ----
  message("\n--- Summary plots (", plot_p_col, ") ---")
  if (do_kegg)
    plot_kegg_summary(out_dir, selected_modules, p_col = plot_p_col,
                      sig_cutoff = plot_sig_cutoff, top_n = plot_top_n,
                      plot_width = plot_width, plot_height = plot_height)
  if (do_go)
    plot_go_summary(out_dir, selected_modules, p_col = plot_p_col,
                    sig_cutoff = plot_sig_cutoff, top_n = plot_top_n,
                    plot_width = plot_width, plot_height = plot_height)
  if (do_reactome)
    plot_reactome_summary(out_dir, selected_modules, p_col = plot_p_col,
                          sig_cutoff = plot_sig_cutoff, top_n = plot_top_n,
                          plot_width = plot_width, plot_height = plot_height)
  if (do_enrichr && length(enrichr_dbs) > 0)
    plot_enrichr_summary(out_dir, selected_modules, enrichr_dbs, p_col = plot_p_col,
                         sig_cutoff = plot_sig_cutoff, top_n = plot_top_n,
                         plot_width = plot_width, plot_height = plot_height)

  # ---- Run parameters ----
  write_lines_safe(
    c(paste0("timestamp\t",              timestamp_now()),
      paste0("project_root\t",           ifelse(is.null(project_root), "", project_root)),
      paste0("variant_label\t",          variant_label),
      paste0("list_stem\t",              list_stem),
      paste0("annotation_dir\t",         annotation_dir),
      paste0("out_dir\t",                out_dir),
      paste0("selection_mode\t",         selection_mode),
      paste0("module_list_file\t",       ifelse(is.null(module_list_file), "", module_list_file)),
      paste0("stats_file\t",             ifelse(is.null(stats_file), "", stats_file)),
      paste0("p_thresh\t",               p_thresh),
      paste0("do_kegg\t",                do_kegg),
      paste0("do_go\t",                  do_go),
      paste0("do_reactome\t",            do_reactome),
      paste0("do_enrichr\t",             do_enrichr),
      paste0("organism\t",               organism),
      paste0("enrichr_dbs\t",            paste(enrichr_dbs, collapse = ",")),
      paste0("use_enrichr_background\t", use_enrichr_background),
      paste0("background_ok\t",          bg_ok),
      paste0("min_bg_genes\t",           min_bg_genes),
      paste0("enrichr_retries\t",        enrichr_retries),
      paste0("enrichr_sleep\t",          enrichr_sleep),
      paste0("plot_p_col\t",             plot_p_col),
      paste0("plot_top_n\t",             plot_top_n),
      paste0("plot_sig_cutoff\t",        plot_sig_cutoff),
      paste0("plot_width\t",             plot_width),
      paste0("plot_height\t",            plot_height),
      paste0("n_selected_modules\t",     length(selected_modules))),
    params_file
  )

  append_log(log_file, "Finished enrichment for variant: ", variant_label,
             if (nzchar(list_stem)) paste0(" / list: ", list_stem) else "")
  invisible(list(selected_modules = selected_modules, out_dir = out_dir))
}

# ============================================================
# 12) Read arguments
# ============================================================
project_root      <- trim_or_null(get_arg("--project_root"))
annotation_dir_v1 <- trim_or_null(get_arg("--annotation_dir_v1"))
annotation_dir_v2 <- trim_or_null(get_arg("--annotation_dir_v2"))
annotation_dir_v3 <- trim_or_null(get_arg("--annotation_dir_v3"))
module_list_input <- trim_or_null(get_arg("--module_list_file"))
stats_file_v1     <- trim_or_null(get_arg("--stats_file_v1"))
stats_file_v2     <- trim_or_null(get_arg("--stats_file_v2"))
stats_file_v3     <- trim_or_null(get_arg("--stats_file_v3"))
trait_list_input  <- trim_or_null(get_arg("--trait_list_file"))
p_thresh          <- as.numeric(get_arg("--p_thresh", "0.05"))

do_kegg     <- to_bool(get_arg("--do_kegg",     "true"),  default = TRUE)
do_go       <- to_bool(get_arg("--do_go",       "true"),  default = TRUE)
do_reactome <- to_bool(get_arg("--do_reactome", "true"),  default = TRUE)
do_enrichr  <- to_bool(get_arg("--do_enrichr",  "false"), default = FALSE)

enrichr_dbs <- split_csv(get_arg(
  "--enrichr_dbs",
  "GO_Biological_Process_2025,GO_Cellular_Component_2025,GO_Molecular_Function_2025"
))

organism               <- trim_or_null(get_arg("--organism", "hsa"))
use_enrichr_background <- to_bool(get_arg("--use_enrichr_background", "true"), TRUE)
min_bg_genes           <- as.integer(get_arg("--min_bg_genes", "200"))
enrichr_retries        <- as.integer(get_arg("--enrichr_retries", "3"))
enrichr_sleep          <- as.numeric(get_arg("--enrichr_sleep", "1"))

plot_p_col      <- trim_or_null(get_arg("--plot_p_col",      "p.adjust"))
plot_top_n      <- as.integer(get_arg( "--plot_top_n",       "10"))
plot_sig_cutoff <- as.numeric(get_arg( "--plot_sig_cutoff",  "0.05"))
plot_width      <- as.numeric(get_arg( "--plot_width",       "14"))
plot_height     <- as.numeric(get_arg( "--plot_height",      "10"))

# ============================================================
# 13) Resolve list inputs into named vectors of file paths
# ============================================================
module_list_files <- resolve_list_input(module_list_input)
trait_list_files  <- resolve_list_input(trait_list_input)

# ============================================================
# 14) Validate top-level arguments
# ============================================================
stop_if_missing(project_root,      "--project_root")
stop_if_missing(annotation_dir_v1, "--annotation_dir_v1")

validate_dir_exists(project_root,      "project_root")
validate_dir_exists(annotation_dir_v1, "annotation_dir_v1")
if (!is.null(annotation_dir_v2)) validate_dir_exists(annotation_dir_v2, "annotation_dir_v2")
if (!is.null(annotation_dir_v3)) validate_dir_exists(annotation_dir_v3, "annotation_dir_v3")

if (!plot_p_col %in% c("p.adjust", "pvalue"))
  stop("--plot_p_col must be 'p.adjust' or 'pvalue'", call. = FALSE)

if (!is.null(module_list_files)) {
  if (!is.null(trait_list_files))
    message("Note: --module_list_file provided; --trait_list_file will be ignored.")
} else {
  if (is.null(trait_list_files))
    stop("Provide --module_list_file OR --trait_list_file with stats_file(s).", call. = FALSE)

  # Trait-list mode requires a matching stats file for every annotation variant
  # that will be run. This avoids vague file.exists(NULL) errors.
  validate_file_exists(stats_file_v1, "stats_file_v1")
  if (!is.null(annotation_dir_v2)) validate_file_exists(stats_file_v2, "stats_file_v2")
  if (!is.null(annotation_dir_v3)) validate_file_exists(stats_file_v3, "stats_file_v3")
}

if (do_enrichr) {
  db_check    <- validate_enrichr_dbs_safe(enrichr_dbs)
  if (length(db_check$bad) > 0)
    warning("Dropping invalid Enrichr DBs: ", paste(db_check$bad, collapse = ", "))
  enrichr_dbs <- db_check$ok
  if (length(enrichr_dbs) == 0)
    stop("After validation, no valid Enrichr DBs remain.", call. = FALSE)
}

setwd(project_root)

cat("project_root: ",      project_root,      "\n", sep = "")
cat("annotation_dir_v1: ", annotation_dir_v1, "\n", sep = "")
if (!is.null(annotation_dir_v2)) cat("annotation_dir_v2: ", annotation_dir_v2, "\n", sep = "")
if (!is.null(annotation_dir_v3)) cat("annotation_dir_v3: ", annotation_dir_v3, "\n", sep = "")
cat("do_kegg: ",         do_kegg,         "\n", sep = "")
cat("do_go: ",           do_go,           "\n", sep = "")
cat("do_reactome: ",     do_reactome,     "\n", sep = "")
cat("do_enrichr: ",      do_enrichr,      "\n", sep = "")
cat("organism: ",        organism,        "\n", sep = "")
cat("plot_p_col: ",      plot_p_col,      "\n", sep = "")
cat("plot_top_n: ",      plot_top_n,      "\n", sep = "")
cat("plot_sig_cutoff: ", plot_sig_cutoff, "\n", sep = "")
if (!is.null(module_list_files)) {
  cat("module list inputs: ", paste(names(module_list_files), basename(module_list_files), sep = "=", collapse = "; "), "\n", sep = "")
}
if (!is.null(trait_list_files)) {
  cat("trait list inputs: ", paste(names(trait_list_files), basename(trait_list_files), sep = "=", collapse = "; "), "\n", sep = "")
}

# ============================================================
# 15) Build the list of (stem, module_list_file, trait_list) combos to iterate
# ============================================================
run_configs <- list()

if (!is.null(module_list_files)) {
  for (i in seq_along(module_list_files)) {
    stem <- names(module_list_files)[i]
    run_configs[[length(run_configs) + 1]] <- list(
      stem             = stem,
      module_list_file = unname(module_list_files[i]),
      trait_list       = NULL
    )
  }
} else {
  for (i in seq_along(trait_list_files)) {
    stem        <- names(trait_list_files)[i]
    trait_list  <- read_trait_list_file(unname(trait_list_files[i]))
    run_configs[[length(run_configs) + 1]] <- list(
      stem             = stem,
      module_list_file = NULL,
      trait_list       = trait_list
    )
  }
}

# ============================================================
# 16) Run all configs × variants
# ============================================================
for (cfg in run_configs) {

  stem_label <- if (nzchar(cfg$stem))
    paste0(" [list: ", cfg$stem, "]") else ""
  message("\n======================================================")
  message("Running config", stem_label)
  message("======================================================")

  common_args <- list(
    module_list_file       = cfg$module_list_file,
    trait_list             = cfg$trait_list,
    p_thresh               = p_thresh,
    do_kegg                = do_kegg,
    do_go                  = do_go,
    do_reactome            = do_reactome,
    do_enrichr             = do_enrichr,
    enrichr_dbs            = enrichr_dbs,
    organism               = organism,
    project_root           = project_root,
    use_enrichr_background = use_enrichr_background,
    min_bg_genes           = min_bg_genes,
    enrichr_retries        = enrichr_retries,
    enrichr_sleep          = enrichr_sleep,
    plot_p_col             = plot_p_col,
    plot_top_n             = plot_top_n,
    plot_sig_cutoff        = plot_sig_cutoff,
    plot_width             = plot_width,
    plot_height            = plot_height,
    list_stem              = cfg$stem
  )

  do.call(run_enrichment_for_variant,
          c(list(annotation_dir = annotation_dir_v1,
                 variant_label  = "v1_all_pcs",
                 stats_file     = stats_file_v1),
            common_args))

  if (!is.null(annotation_dir_v2)) {
    do.call(run_enrichment_for_variant,
            c(list(annotation_dir = annotation_dir_v2,
                   variant_label  = "v2_exclude_protected_pcs",
                   stats_file     = stats_file_v2),
              common_args))
  }

  if (!is.null(annotation_dir_v3)) {
    do.call(run_enrichment_for_variant,
            c(list(annotation_dir = annotation_dir_v3,
                   variant_label  = "v3_technical_pcs_only",
                   stats_file     = stats_file_v3),
              common_args))
  }
}

cat("\nScript 12b complete: Enrichment finished\n")