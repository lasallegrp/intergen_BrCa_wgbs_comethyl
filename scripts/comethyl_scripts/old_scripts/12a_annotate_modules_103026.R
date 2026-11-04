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
#     --annotation_mode auto \
#     --cpg_island_file /path/to/local/cpgIslandExt.txt.gz \
#     [--add_genomic_location TRUE] \
#     [--helper_file /path/to/helper.R]   # optional; defaults to
#                                          # helper.R alongside this script
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
#     Module_Genic_Enrichment.tsv / .xlsx / .pdf   (module x genic-category Fisher's tests + heatmap)
#     Module_CpG_Enrichment.tsv / .xlsx / .pdf     (module x CpG-category Fisher's tests + heatmap)
#     session_info.txt
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
# 1) Shared helpers (CLI parsing, logging, validation, annotation core)
# ============================================================
.get_arg_bootstrap <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) NULL
)
if (is.null(script_dir) || !nzchar(script_dir)) script_dir <- getwd()

# Default: helper.R lives alongside this script. Override with
# --helper_file /full/path/to/helper.R if it lives elsewhere.
helper_file <- .get_arg_bootstrap("--helper_file", file.path(script_dir, "helper.R"))
if (!file.exists(helper_file)) {
  stop("helper.R not found at: ", helper_file,
       "\nPass --helper_file /path/to/helper.R if it lives elsewhere.", call. = FALSE)
}
source(helper_file)

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
# 6) Module/gene summary helpers
# ============================================================
# (offline_nearest_gene, annotate_offline_only, annotate_regions_safe
#  now come from helper.R)
# ============================================================
# extract_regions_from_module_object() now comes from helper.R

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
                                       sample_info_path,
                                       sample_id_col,
                                       genome,
                                       annotation_mode,
                                       project_root,
                                       add_genomic_location = TRUE,
                                       cpg_island_file = NULL,
                                       txdb_pkg = NULL) {

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
  append_log(log_file, "add_genomic_location: ", add_genomic_location)
  append_log(log_file, "txdb_pkg: ", ifelse(is.null(txdb_pkg), paste0("(auto-resolved from genome=", genome, ")"), txdb_pkg))
  append_log(log_file, "cpg_island_file: ", ifelse(is.null(cpg_island_file), "(not provided)", cpg_island_file))

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
      verbose = TRUE,
      add_genomic_location = add_genomic_location,
      cpg_island_file = cpg_island_file,
      txdb_pkg = txdb_pkg
    )
  )

  append_log(log_file, "Annotated regions rows: ", nrow(annotated_regions))
  if ("genomic_location" %in% colnames(annotated_regions)) {
    append_log(log_file, "Regions with genomic_location: ",
               sum(!is.na(annotated_regions$genomic_location)))
  }
  if ("cpg_context" %in% colnames(annotated_regions)) {
    append_log(log_file, "Regions with cpg_context: ",
               sum(!is.na(annotated_regions$cpg_context)))
  }

  openxlsx::write.xlsx(annotated_regions, annotated_xlsx, rowNames = FALSE)

  # --------------------------------------------------------
  # Module-level annotation enrichment (genic + CpG), only when the
  # required columns are present -- i.e. offline path ran with
  # add_genomic_location/cpg_island_file. Skipped silently (with a
  # log note) for GREAT-based output, which doesn't have these columns.
  # --------------------------------------------------------
  if ("genomic_location" %in% colnames(annotated_regions)) {
    genic_enrichment <- moduleGenicEnrichment(annotated_regions, verbose = TRUE)
    write.table(genic_enrichment, file.path(out_dir, "Module_Genic_Enrichment.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    openxlsx::write.xlsx(genic_enrichment, file.path(out_dir, "Module_Genic_Enrichment.xlsx"), rowNames = FALSE)
    plotModuleAnnotationEnrichment(
      genic_enrichment,
      title = paste0("Genic Region Enrichment by Module (", variant_label, ")"),
      file = file.path(out_dir, "Module_Genic_Enrichment.pdf")
    )
    append_log(log_file, "Wrote genic enrichment table + plot (",
               nrow(genic_enrichment), " module x category tests)")
  } else {
    append_log(log_file, "Skipping genic enrichment (no genomic_location column)")
  }

  if (all(c("CpG.Island", "CpG.Shore", "CpG.Shelf", "Open.Sea") %in% colnames(annotated_regions))) {
    cpg_enrichment <- moduleCpGEnrichment(annotated_regions, verbose = TRUE)
    write.table(cpg_enrichment, file.path(out_dir, "Module_CpG_Enrichment.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    openxlsx::write.xlsx(cpg_enrichment, file.path(out_dir, "Module_CpG_Enrichment.xlsx"), rowNames = FALSE)
    plotModuleAnnotationEnrichment(
      cpg_enrichment,
      title = paste0("CpG Context Enrichment by Module (", variant_label, ")"),
      categoryOrder = c("Island", "Shore", "Shelf", "OpenSea"),
      file = file.path(out_dir, "Module_CpG_Enrichment.pdf")
    )
    append_log(log_file, "Wrote CpG enrichment table + plot (",
               nrow(cpg_enrichment), " module x category tests)")
  } else {
    append_log(log_file, "Skipping CpG enrichment (no CpG.Island/Shore/Shelf/Open.Sea columns)")
  }

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
    paste0("add_genomic_location\t", add_genomic_location),
    paste0("cpg_island_file\t", ifelse(is.null(cpg_island_file), "", cpg_island_file)),
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

  session_info_file <- file.path(out_dir, "session_info.txt")
  write_session_info(
    session_info_file,
    extra_files = list(
      modules_path = modules_path,
      sample_info = sample_info_path,
      cpg_island_file = cpg_island_file
    )
  )
  append_log(log_file, "Wrote: ", session_info_file)

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
cpg_island_file <- trim_or_null(get_arg("--cpg_island_file"))
add_genomic_location <- parse_bool(get_arg("--add_genomic_location", "TRUE"), "--add_genomic_location")
txdb_pkg        <- trim_or_null(get_arg("--txdb_pkg"))  # override; default NULL = auto-resolve from --genome

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
if (!is.null(cpg_island_file)) validate_file_exists(cpg_island_file, "cpg_island_file")

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
  sample_info_path = sample_info,
  sample_id_col = sample_id_col,
  genome = genome,
  annotation_mode = annotation_mode,
  project_root = project_root,
  add_genomic_location = add_genomic_location,
  cpg_island_file = cpg_island_file,
  txdb_pkg = txdb_pkg
)

if (!is.null(modules_v2)) {
  run_annotation_for_variant(
    modules_path = modules_v2,
    variant_label = "v2_exclude_protected_pcs",
    sample_info_df = sample_info_df,
    sample_info_path = sample_info,
    sample_id_col = sample_id_col,
    genome = genome,
    annotation_mode = annotation_mode,
    project_root = project_root,
    add_genomic_location = add_genomic_location,
    cpg_island_file = cpg_island_file,
    txdb_pkg = txdb_pkg
  )
}

if (!is.null(modules_v3)) {
  run_annotation_for_variant(
    modules_path = modules_v3,
    variant_label = "v3_technical_pcs_only",
    sample_info_df = sample_info_df,
    sample_info_path = sample_info,
    sample_id_col = sample_id_col,
    genome = genome,
    annotation_mode = annotation_mode,
    project_root = project_root,
    add_genomic_location = add_genomic_location,
    cpg_island_file = cpg_island_file,
    txdb_pkg = txdb_pkg
  )
}

message("Script 12a complete: Annotate Modules finished")