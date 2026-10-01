#!/usr/bin/env Rscript

# ============================================================
# 12A_annotate_modules.R
#
# Purpose:
#   Reproducible annotation step for comethyl module objects.
#
# What this script does:
#   - loads one or more module objects
#   - extracts module regions
#   - annotates all regions (all modules, including grey)
#   - writes annotated region tables
#   - writes per-module gene lists
#   - writes module-level summary tables
#   - writes background/testable gene universe
#   - writes run logs
#
# What this script does NOT do:
#   - module significance filtering
#   - KEGG enrichment
#   - Enrichr enrichment
#
# Expected module object:
#   A comethyl module RDS with at least:
#     $regions
#   where regions contains columns:
#     RegionID, chr, start, end, module
#
# Recommended usage:
#   Rscript 12A_annotate_modules.R \
#     --project_root /path/to/project \
#     --modules_v1 /path/to/v1/modules/Modules.rds \
#     --modules_v2 /path/to/v2/modules/Modules.rds \
#     --modules_v3 /path/to/v3/modules/Modules.rds \
#     --sample_info /path/to/sample_info.xlsx \
#     --genome hg38 \
#     --annotation_mode auto
#
# Output structure per variant:
#   <project_root>/comethyl_output/12_annotation/<cpg_label>/<region_label>/<variant_name>/
#     Annotated_Regions.tsv
#     Annotated_Regions.xlsx
#     Module_Gene_List.tsv
#     Module_Gene_Summary.tsv
#     Module_Gene_Summary.xlsx
#     Background_Genes.txt
#     Background_Genes.tsv
#     run_parameters.txt
#     run_log.txt
# ============================================================

message("Starting Script 12a")

suppressPackageStartupMessages({
  library(comethyl)
  library(dplyr)
  library(openxlsx)
  library(readxl)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(AnnotationDbi)
  library(stringr)
})

# ============================================================
# 1) Base argument parsing
# ============================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

trim_or_null <- function(x) {
  if (is.null(x) || is.na(x)) return(NULL)
  x <- trimws(x)
  if (!nzchar(x)) return(NULL)
  x
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

append_log <- function(logfile, ...) {
  txt <- paste0("[", timestamp_now(), "] ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = logfile, append = TRUE)
}

write_lines_safe <- function(x, file) {
  writeLines(as.character(x), con = file, useBytes = TRUE)
}

# ============================================================
# 2) Validation helpers
# ============================================================
stop_if_missing <- function(x, label) {
  if (is.null(x) || !nzchar(x)) stop("Missing required argument: ", label, call. = FALSE)
}

validate_file_exists <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

validate_regions_df <- function(regions, source_label = "regions") {
  req <- c("RegionID", "chr", "start", "end", "module")
  missing_cols <- setdiff(req, colnames(regions))
  if (length(missing_cols) > 0) {
    stop(
      source_label, " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

# ============================================================
# 3) Sample info loader
# ============================================================
load_sample_info <- function(sample_info, sample_id_col = NULL) {
  ext <- tolower(tools::file_ext(sample_info))

  if (ext %in% c("xlsx", "xls")) {
    df <- openxlsx::read.xlsx(sample_info, rowNames = FALSE)
  } else if (ext %in% c("csv")) {
    df <- read.csv(sample_info, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    df <- read.delim(sample_info, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported sample_info format: ", sample_info, call. = FALSE)
  }

  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(df)) {
      stop("sample_id_col not found in sample_info: ", sample_id_col, call. = FALSE)
    }
  }

  df
}

# ============================================================
# 4) Path derivation helper
# ============================================================
derive_pipeline_dirs_from_modules <- function(modules_rds, project_root, step_name = "12_annotation") {
  variant_dir  <- dirname(modules_rds)
  region_dir   <- dirname(variant_dir)
  variant_name <- basename(variant_dir)
  region_label <- basename(region_dir)
  cpg_label    <- basename(dirname(region_dir))

  pipeline_root <- file.path(project_root, "comethyl_output")
  step_dir <- file.path(pipeline_root, step_name)
  out_dir <- file.path(step_dir, cpg_label, region_label, variant_name)

  safe_dir_create(out_dir)

  list(
    pipeline_root = pipeline_root,
    step_dir = step_dir,
    out_dir = out_dir,
    variant_name = variant_name,
    region_label = region_label,
    cpg_label = cpg_label
  )
}

# ============================================================
# 5) Annotation helpers
# ============================================================
offline_nearest_gene <- function(gr, verbose = TRUE) {
  have_ensdb <- requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)
  have_txdb  <- requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)
  have_org   <- requireNamespace("org.Hs.eg.db", quietly = TRUE)

  if (!have_org) {
    stop("org.Hs.eg.db is required for offline annotation.", call. = FALSE)
  }

  if (have_ensdb) {
    if (verbose) message("[offline] Using EnsDb.Hsapiens.v86")
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    seqs <- unique(as.character(GenomeInfoDb::seqnames(gr)))
    seqs <- seqs[seqs %in% GenomeInfoDb::seqlevels(edb)]

    genes <- ensembldb::genes(edb, filter = AnnotationFilter::SeqNameFilter(seqs))
    ggr <- GenomicRanges::GRanges(genes)

    ens_ids <- genes$gene_id
    map <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = ens_ids,
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "ENTREZID")
    )

    ggr$SYMBOL   <- map$SYMBOL[match(genes$gene_id, map$ENSEMBL)]
    ggr$ENTREZID <- map$ENTREZID[match(genes$gene_id, map$ENSEMBL)]

  } else if (have_txdb) {
    if (verbose) message("[offline] Using TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    ggr  <- GenomicFeatures::genes(txdb)

    map <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = ggr$gene_id,
      keytype = "ENTREZID",
      columns = c("SYMBOL")
    )

    ggr$ENTREZID <- ggr$gene_id
    ggr$SYMBOL   <- map$SYMBOL[match(ggr$ENTREZID, map$ENTREZID)]

  } else {
    stop(
      "Install one of: EnsDb.Hsapiens.v86 OR TxDb.Hsapiens.UCSC.hg38.knownGene",
      call. = FALSE
    )
  }

  hit <- GenomicRanges::distanceToNearest(gr, ggr, ignore.strand = TRUE)
  ng  <- ggr[S4Vectors::subjectHits(hit)]

  data.frame(
    chr = as.character(GenomeInfoDb::seqnames(gr))[S4Vectors::queryHits(hit)],
    start = as.integer(S4Vectors::start(gr))[S4Vectors::queryHits(hit)],
    end = as.integer(S4Vectors::end(gr))[S4Vectors::queryHits(hit)],
    gene_symbol = as.character(ng$SYMBOL),
    gene_entrezID = as.character(ng$ENTREZID),
    stringsAsFactors = FALSE
  )
}

annotate_offline_only <- function(regions_df, genome = "hg38", file_txt = NULL, verbose = TRUE) {
  validate_regions_df(regions_df, "regions_df")

  gr <- GenomicRanges::GRanges(
    seqnames = regions_df$chr,
    ranges   = IRanges::IRanges(start = regions_df$start, end = regions_df$end),
    RegionID = regions_df$RegionID
  )

  ng <- offline_nearest_gene(gr, verbose = verbose)

  out <- regions_df %>%
    left_join(
      ng %>% dplyr::select(chr, start, end, gene_symbol, gene_entrezID),
      by = c("chr", "start", "end")
    ) %>%
    mutate(
      gene_description = NA_character_,
      gene_ensemblID   = NA_character_
    )

  if (!is.null(file_txt) && nzchar(file_txt)) {
    write.table(out, file = file_txt, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  out
}

annotate_regions_safe <- function(regions_df,
                                  genome = "hg38",
                                  annotation_mode = c("auto", "great", "offline"),
                                  file_txt = NULL,
                                  verbose = TRUE) {
  annotation_mode <- match.arg(annotation_mode)

  if (annotation_mode == "offline") {
    return(annotate_offline_only(regions_df, genome = genome, file_txt = file_txt, verbose = verbose))
  }

  if (annotation_mode == "great") {
    if (verbose) message("[annotate] Using comethyl::annotateModule() only")
    return(comethyl::annotateModule(regions_df, genome = genome, file = file_txt))
  }

  tryCatch(
    {
      if (verbose) message("[annotate] Trying comethyl::annotateModule()")
      comethyl::annotateModule(regions_df, genome = genome, file = file_txt)
    },
    error = function(e) {
      message("[annotate] GREAT-based annotation failed: ", conditionMessage(e))
      message("[annotate] Falling back to offline nearest-gene annotation.")
      annotate_offline_only(regions_df, genome = genome, file_txt = file_txt, verbose = verbose)
    }
  )
}

# ============================================================
# 6) Module/gene summary helpers
# ============================================================
extract_regions_from_module_object <- function(obj, label = "module object") {
  if (is.list(obj) && "regions" %in% names(obj)) {
    regions <- obj$regions
  } else if (is.data.frame(obj)) {
    regions <- obj
  } else {
    stop(
      "Could not extract regions from ", label,
      ". Expected either a list with $regions or a data.frame.",
      call. = FALSE
    )
  }

  validate_regions_df(regions, label)
  regions
}

make_module_gene_list_table <- function(annotated_regions) {
  validate_regions_df(annotated_regions, "annotated_regions")

  if (!"gene_symbol" %in% colnames(annotated_regions)) {
    stop("annotated_regions must contain gene_symbol column.", call. = FALSE)
  }

  annotated_regions %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::select(module, gene_symbol) %>%
    dplyr::distinct() %>%
    dplyr::arrange(module, gene_symbol)
}

make_module_summary_table <- function(annotated_regions) {
  validate_regions_df(annotated_regions, "annotated_regions")

  gene_col_present <- "gene_symbol" %in% colnames(annotated_regions)

  annotated_regions %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(
      n_regions = dplyr::n(),
      n_unique_genes = if (gene_col_present) length(unique(na.omit(gene_symbol[gene_symbol != ""]))) else NA_integer_,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_regions), module)
}

make_background_gene_table <- function(annotated_regions) {
  if (!"gene_symbol" %in% colnames(annotated_regions)) {
    stop("annotated_regions must contain gene_symbol column.", call. = FALSE)
  }

  tibble::tibble(
    gene_symbol = sort(unique(na.omit(as.character(annotated_regions$gene_symbol))))
  ) %>%
    dplyr::filter(gene_symbol != "")
}

# ============================================================
# 7) Per-variant runner
# ============================================================
run_annotation_for_variant <- function(modules_path,
                                       variant_label,
                                       sample_info_df,
                                       sample_id_col,
                                       genome,
                                       annotation_mode,
                                       project_root) {

  validate_file_exists(modules_path, paste0("modules_", variant_label))

  dir_info <- derive_pipeline_dirs_from_modules(
    modules_rds = modules_path,
    project_root = project_root,
    step_name = "12_annotation"
  )

  out_dir <- dir_info$out_dir

  log_file <- file.path(out_dir, "run_log.txt")
  params_file <- file.path(out_dir, "run_parameters.txt")

  append_log(log_file, "Starting annotation for variant: ", variant_label)
  append_log(log_file, "modules_path: ", modules_path)
  append_log(log_file, "pipeline_root: ", dir_info$pipeline_root)
  append_log(log_file, "step_dir: ", dir_info$step_dir)
  append_log(log_file, "cpg_label: ", dir_info$cpg_label)
  append_log(log_file, "region_label: ", dir_info$region_label)
  append_log(log_file, "variant_name: ", dir_info$variant_name)
  append_log(log_file, "out_dir: ", out_dir)
  append_log(log_file, "genome: ", genome)
  append_log(log_file, "annotation_mode: ", annotation_mode)

  if (!is.null(sample_info_df)) {
    append_log(log_file, "sample_info rows: ", nrow(sample_info_df), "; cols: ", ncol(sample_info_df))
    if (!is.null(sample_id_col)) {
      append_log(log_file, "sample_id_col: ", sample_id_col)
    } else {
      append_log(log_file, "sample_id_col: (not provided)")
    }
  }

  obj <- readRDS(modules_path)
  regions <- extract_regions_from_module_object(obj, label = modules_path)

  append_log(log_file, "Loaded regions: ", nrow(regions))
  append_log(log_file, "Unique modules: ", length(unique(as.character(regions$module))))

  annotated_tsv  <- file.path(out_dir, "Annotated_Regions.tsv")
  annotated_xlsx <- file.path(out_dir, "Annotated_Regions.xlsx")
  gene_list_tsv  <- file.path(out_dir, "Module_Gene_List.tsv")
  gene_sum_tsv   <- file.path(out_dir, "Module_Gene_Summary.tsv")
  gene_sum_xlsx  <- file.path(out_dir, "Module_Gene_Summary.xlsx")
  bg_txt         <- file.path(out_dir, "Background_Genes.txt")
  bg_tsv         <- file.path(out_dir, "Background_Genes.tsv")

  annotated_regions <- suppressWarnings(
    annotate_regions_safe(
      regions_df = regions,
      genome = genome,
      annotation_mode = annotation_mode,
      file_txt = annotated_tsv,
      verbose = TRUE
    )
  )

  append_log(log_file, "Annotated regions rows: ", nrow(annotated_regions))

  openxlsx::write.xlsx(annotated_regions, annotated_xlsx, rowNames = FALSE)

  module_gene_list <- make_module_gene_list_table(annotated_regions)
  write.table(module_gene_list, gene_list_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  module_summary <- make_module_summary_table(annotated_regions)
  write.table(module_summary, gene_sum_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  openxlsx::write.xlsx(module_summary, gene_sum_xlsx, rowNames = FALSE)

  background_genes <- make_background_gene_table(annotated_regions)
  write_lines_safe(background_genes$gene_symbol, bg_txt)
  write.table(background_genes, bg_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  append_log(log_file, "Module gene list rows: ", nrow(module_gene_list))
  append_log(log_file, "Module summary rows: ", nrow(module_summary))
  append_log(log_file, "Background genes: ", nrow(background_genes))

  params <- c(
    paste0("timestamp\t", timestamp_now()),
    paste0("project_root\t", project_root),
    paste0("modules_path\t", modules_path),
    paste0("genome\t", genome),
    paste0("annotation_mode\t", annotation_mode),
    paste0("sample_info_provided\t", !is.null(sample_info_df)),
    paste0("sample_id_col\t", ifelse(is.null(sample_id_col), "", sample_id_col)),
    paste0("pipeline_root\t", dir_info$pipeline_root),
    paste0("step_dir\t", dir_info$step_dir),
    paste0("cpg_label\t", dir_info$cpg_label),
    paste0("region_label\t", dir_info$region_label),
    paste0("variant_name\t", dir_info$variant_name),
    paste0("out_dir\t", out_dir),
    paste0("n_regions\t", nrow(regions)),
    paste0("n_modules\t", length(unique(as.character(regions$module)))),
    paste0("n_annotated_regions\t", nrow(annotated_regions)),
    paste0("n_background_genes\t", nrow(background_genes))
  )
  write_lines_safe(params, params_file)

  append_log(log_file, "Finished annotation for variant: ", variant_label)

  invisible(list(
    annotated_regions = annotated_regions,
    module_gene_list = module_gene_list,
    module_summary = module_summary,
    background_genes = background_genes,
    out_dir = out_dir
  ))
}

# ============================================================
# 8) Read arguments
# ============================================================
project_root    <- trim_or_null(get_arg("--project_root"))
sample_info     <- trim_or_null(get_arg("--sample_info"))
sample_id_col   <- trim_or_null(get_arg("--sample_id_col"))

modules_v1      <- trim_or_null(get_arg("--modules_v1"))
modules_v2      <- trim_or_null(get_arg("--modules_v2"))
modules_v3      <- trim_or_null(get_arg("--modules_v3"))

genome          <- trim_or_null(get_arg("--genome", "hg38"))
annotation_mode <- trim_or_null(get_arg("--annotation_mode", "auto"))

# ============================================================
# 9) Validate required inputs
# ============================================================
stop_if_missing(project_root, "--project_root")
stop_if_missing(sample_info, "--sample_info")
stop_if_missing(modules_v1, "--modules_v1")

if (!dir.exists(project_root)) {
  stop("project_root not found: ", project_root, call. = FALSE)
}
validate_file_exists(sample_info, "sample_info")
validate_file_exists(modules_v1, "modules_v1")

if (!is.null(modules_v2)) validate_file_exists(modules_v2, "modules_v2")
if (!is.null(modules_v3)) validate_file_exists(modules_v3, "modules_v3")

annotation_mode <- match.arg(annotation_mode, choices = c("auto", "great", "offline"))

setwd(project_root)

AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(project_root, ".cache")
)

# ============================================================
# 10) Load sample info once
# ============================================================
sample_info_df <- load_sample_info(sample_info, sample_id_col = sample_id_col)

cat("Loaded sample_info: ", sample_info, "\n", sep = "")
cat("sample_info dimensions: ", nrow(sample_info_df), " x ", ncol(sample_info_df), "\n", sep = "")

if (!is.null(sample_id_col)) {
  cat("Using sample_id_col: ", sample_id_col, "\n", sep = "")
} else {
  cat("No sample_id_col provided.\n")
}

# ============================================================
# 11) Run variants
# ============================================================
run_annotation_for_variant(
  modules_path = modules_v1,
  variant_label = "v1_all_pcs",
  sample_info_df = sample_info_df,
  sample_id_col = sample_id_col,
  genome = genome,
  annotation_mode = annotation_mode,
  project_root = project_root
)

if (!is.null(modules_v2)) {
  run_annotation_for_variant(
    modules_path = modules_v2,
    variant_label = "v2_exclude_protected_pcs",
    sample_info_df = sample_info_df,
    sample_id_col = sample_id_col,
    genome = genome,
    annotation_mode = annotation_mode,
    project_root = project_root
  )
}

if (!is.null(modules_v3)) {
  run_annotation_for_variant(
    modules_path = modules_v3,
    variant_label = "v3_technical_pcs_only",
    sample_info_df = sample_info_df,
    sample_id_col = sample_id_col,
    genome = genome,
    annotation_mode = annotation_mode,
    project_root = project_root
  )
}

message("Script 12a complete: Annotate Modules finished")