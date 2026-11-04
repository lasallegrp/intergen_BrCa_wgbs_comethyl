#!/usr/bin/env Rscript
# ================================================================
# SCRIPT 04: PC Derivation and PC-Trait Diagnostics
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   - Load region-level methylation from Script 03
#   - Load analysis-ready sample information
#   - Align samples between methylation and sample info
#   - Derive top PCs/SVs using sva::num.sv + SVD
#   - Allow explicit selection of the model used by sva::num.sv:
#       * intercept_only  = original comethyl-style PC/SV derivation
#       * protected       = protected-trait-aware PC/SV derivation
#   - Compute PC-trait association statistics using Pearson and bicor
#   - Generate PC-trait diagnostic heatmaps and trait dendrogram
#   - Save reusable PC/SV outputs for downstream methylation adjustment
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --input_meth   : path to Region_Methylation.rds from Script 03
#   --sample_info  : path to sample information file (.xlsx, .csv, .tsv, .txt)
#
# OPTIONAL INPUTS
#   --sample_id_col         : column in sample_info containing sample IDs
#                             [default = use existing rownames]
#
#   --protected_traits_file : text file with one protected trait per line.
#                             These traits are used for:
#                               1. selected-trait diagnostic heatmaps
#                               2. PC-trait association summaries
#                               3. protected-model PC/SV derivation when
#                                  --pc_model protected is selected
#
#   --technical_traits_file : text file with one technical trait per line.
#                             These traits are used for diagnostics and optional
#                             selected-trait plots, but are not used in the
#                             protected model for num.sv().
#
#   --pc_model              : model used by sva::num.sv()
#                             Allowed values:
#                               intercept_only
#                               protected
#
#                             [default = intercept_only]
#
#                             intercept_only:
#                               Uses mod = ~1 for sva::num.sv().
#                               This preserves the original comethyl workflow:
#                               derive SVs using an intercept-only model, then
#                               downstream v1 adjusts methylation for all
#                               selected PCs/SVs.
#
#                             protected:
#                               Uses the traits in --protected_traits_file as
#                               the model matrix for sva::num.sv().
#                               This is intended for protected-trait-aware v2
#                               analyses where exposure/outcome-associated PCs
#                               are identified and excluded from adjustment.
#
#   --threads               : number of WGCNA threads [default = 4]
#   --top_trait_n           : number of top traits for overall heatmap [default = 50]
#   --top_cell_n            : number of top PC-trait cells for association heatmap [default = 250]
#   --primary_cor_method    : primary correlation method for diagnostic plots:
#                             bicor or pearson [default = bicor]
#   --selected_trait_mode   : protected, protected_and_technical, none [default = auto]
#   --max_pcs_plot          : number of PCs to show in variance explained plot [default = 20]
#   --plot_width            : plot width in inches [default = 12]
#   --plot_height           : plot height in inches [default = 10]
#   --plot_dpi              : plot DPI [default = 600]
#   --base_size             : base text size for plots [default = 11]
#   --axis_text_size        : axis text size for plots [default = 8]
#
# OUTPUTS
#   project_root/comethyl_output/04_pc_diagnostics/<cpg_label>/<region_label>/<pc_model_label>/
#
#   where <pc_model_label> is:
#       intercept_only
#       protected_model
#
#   Files written:
#       PCs.rds
#           - matrix of selected PCs/SVs
#           - rows = samples
#           - columns = PC_1 ... PC_N selected by sva::num.sv()
#
#       PC_Metadata.tsv
#           - singular values, variance explained, cumulative variance explained
#           - flags which PCs are selected as surrogate variables
#
#       Variance_Explained_Top20_PCs.pdf
#           - PCs selected as surrogate variables by num.sv are shown in red
#           - remaining PCs shown in grey
#
#       PC_Trait_Correlation_Stats_Pearson.tsv
#       PC_Trait_Correlation_Stats_Bicor.tsv
#           - association statistics between selected PCs/SVs and numeric traits
#
#       PC_Trait_Heatmap_SelectedTraits.pdf
#           - generated if selected traits are available
#           - x-axis labels use "PC_N [SV]" for num.sv-selected PCs
#
#       PC_Trait_Heatmap_Top250Associations.pdf
#       PC_Trait_Heatmap_Top50Traits_Overall.pdf
#       Trait_Dendrogram_From_PC_Associations.pdf
#
#       sample_info_non_numeric_columns_removed.txt
#       sample_info_all_na_numeric_columns_removed.txt
#       sample_info_zero_variance_numeric_columns_removed.txt
#
#       traits_requested_protected.txt
#       traits_found_protected.txt
#       traits_missing_protected.txt
#       traits_requested_technical.txt
#       traits_found_technical.txt
#       traits_missing_technical.txt
#
#       sample_alignment/samples_only_in_methylation.txt
#       sample_alignment/samples_only_in_sample_info.txt
#
#       run_parameters.txt
#           - records pc_model, pc_model_label, model formula, sample counts,
#             number of selected PCs/SVs, trait counts, and run date
#
# NOTES
#   - sample_info is assumed to be analysis-ready before entering the pipeline.
#   - No project-specific trait recoding is performed here.
#   - The default --pc_model is intercept_only to preserve the original
#     comethyl behavior.
#   - v1 should use PCs derived with --pc_model intercept_only and should adjust
#     methylation for all selected PCs/SVs.
#   - v2 should use PCs derived with --pc_model protected, then exclude PCs
#     associated with protected exposure/outcome traits before methylation
#     adjustment.
#   - Changing --protected_traits_file can change the protected-model PCs/SVs,
#     but should not change intercept-only PC/SV derivation.
#   - Protected traits are still useful in intercept_only mode for diagnostic
#     PC-trait correlations and selected-trait heatmaps.
#   - TERMINOLOGY: All components from SVD are called PCs. The subset identified
#     by sva::num.sv() as statistically meaningful are labelled "SV"
#     (surrogate variable) in plots to distinguish them from the full PC set.
#     The [SV] tag appears on x-axis labels in heatmaps and SVs are highlighted
#     in red in the variance explained plot.
# ================================================================
# ================================================================
# SCRIPT 04: PC Derivation and PC-Trait Diagnostics
#
# Key update:
#   --pc_model intercept_only  = original comethyl behavior, mod = ~1 [default]
#   --pc_model protected       = use protected traits in sva::num.sv model
#
# Outputs are separated by PC model:
#   comethyl_output/04_pc_diagnostics/<cpg_label>/<region_label>/intercept_only/
#   comethyl_output/04_pc_diagnostics/<cpg_label>/<region_label>/protected_model/
# ================================================================

message("Starting Script 04")

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(comethyl)
  library(sva)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(AnnotationHub)
})

# ------------------------------------------------------------
# Load helper functions
# ------------------------------------------------------------
script_file_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_file_arg) == 0) stop("Could not determine script path from commandArgs().")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found next to this script: ", helper_file)
source(helper_file)

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),
  make_option("--input_meth", type = "character",
              help = "Path to Region_Methylation.rds from Script 03"),
  make_option("--sample_info", type = "character",
              help = "Path to sample information file (.xlsx, .csv, .tsv, .txt)"),
  make_option("--sample_id_col", type = "character", default = NULL,
              help = "Column in sample_info containing sample IDs [default = rownames]"),
  make_option("--protected_traits_file", type = "character", default = NULL,
              help = "Optional text file with one protected trait per line"),
  make_option("--technical_traits_file", type = "character", default = NULL,
              help = "Optional text file with one technical trait per line"),
  make_option("--pc_model", type = "character", default = "intercept_only",
              help = "Model used by sva::num.sv(): intercept_only or protected [default = intercept_only]"),
  make_option("--threads", type = "integer", default = 4,
              help = "Number of WGCNA threads [default = 4]"),
  make_option("--top_trait_n", type = "integer", default = 50,
              help = "Number of top traits for overall heatmap [default = 50]"),
  make_option("--top_cell_n", type = "integer", default = 250,
              help = "Number of top PC-trait cells for association heatmap [default = 250]"),
  make_option("--primary_cor_method", type = "character", default = "bicor",
              help = "Primary correlation method for diagnostics: bicor or pearson [default = bicor]"),
  make_option("--selected_trait_mode", type = "character", default = NULL,
              help = "Trait mode for selected-trait heatmap: protected, protected_and_technical, or none [default = auto]"),
  make_option("--max_pcs_plot", type = "integer", default = 20,
              help = "Number of PCs to display in variance explained plot [default = 20]"),
  make_option("--plot_width", type = "double", default = 12,
              help = "Plot width in inches [default = 12]"),
  make_option("--plot_height", type = "double", default = 10,
              help = "Plot height in inches [default = 10]"),
  make_option("--plot_dpi", type = "integer", default = 600,
              help = "Plot DPI [default = 600]"),
  make_option("--base_size", type = "double", default = 11,
              help = "Base plot text size [default = 11]"),
  make_option("--axis_text_size", type = "double", default = 8,
              help = "Axis text size [default = 8]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$input_meth))   stop("--input_meth is required")
if (is.null(opt$sample_info))  stop("--sample_info is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$input_meth))  stop("Input methylation file not found: ", opt$input_meth)
if (!file.exists(opt$sample_info)) stop("Sample info file not found: ", opt$sample_info)
if (!is.null(opt$protected_traits_file) && !file.exists(opt$protected_traits_file)) {
  stop("Protected traits file not found: ", opt$protected_traits_file)
}
if (!is.null(opt$technical_traits_file) && !file.exists(opt$technical_traits_file)) {
  stop("Technical traits file not found: ", opt$technical_traits_file)
}
if (opt$threads < 1)      stop("--threads must be >= 1")
if (opt$top_trait_n < 1)  stop("--top_trait_n must be >= 1")
if (opt$top_cell_n < 1)   stop("--top_cell_n must be >= 1")
if (opt$max_pcs_plot < 1) stop("--max_pcs_plot must be >= 1")

opt$primary_cor_method <- tolower(opt$primary_cor_method)
if (!opt$primary_cor_method %in% c("bicor", "pearson")) {
  stop("--primary_cor_method must be 'bicor' or 'pearson'")
}

opt$pc_model <- tolower(opt$pc_model)
if (!opt$pc_model %in% c("intercept_only", "protected")) {
  stop("--pc_model must be either 'intercept_only' or 'protected'")
}

if (!is.null(opt$selected_trait_mode) &&
    !tolower(opt$selected_trait_mode) %in% c("protected", "protected_and_technical", "none")) {
  stop("--selected_trait_mode must be one of: protected, protected_and_technical, none")
}

# ------------------------------------------------------------
# Configure cache and threads
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption("CACHE", value = file.path(opt$project_root, ".cache"))
options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads(nThreads = opt$threads)

# ------------------------------------------------------------
# Derive output directory from methylation lineage
# ------------------------------------------------------------
input_dir    <- dirname(opt$input_meth)
region_label <- basename(input_dir)
cpg_label    <- basename(dirname(input_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir      <- file.path(pipeline_root, "04_pc_diagnostics")
pc_model_label <- ifelse(opt$pc_model == "intercept_only", "intercept_only", "protected_model")
out_dir <- file.path(step_dir, cpg_label, region_label, pc_model_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)
message("PC model requested: ", opt$pc_model)

# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------
pcs_rds               <- file.path(out_dir, "PCs.rds")
pc_meta_tsv           <- file.path(out_dir, "PC_Metadata.tsv")
var_plot_pdf          <- file.path(out_dir, paste0("Variance_Explained_Top", opt$max_pcs_plot, "_PCs.pdf"))
pc_trait_pearson_tsv  <- file.path(out_dir, "PC_Trait_Correlation_Stats_Pearson.tsv")
pc_trait_bicor_tsv    <- file.path(out_dir, "PC_Trait_Correlation_Stats_Bicor.tsv")
heatmap_selected_pdf  <- file.path(out_dir, "PC_Trait_Heatmap_SelectedTraits.pdf")
heatmap_topcells_pdf  <- file.path(out_dir, paste0("PC_Trait_Heatmap_Top", opt$top_cell_n, "Associations.pdf"))
heatmap_toptraits_pdf <- file.path(out_dir, paste0("PC_Trait_Heatmap_Top", opt$top_trait_n, "Traits_Overall.pdf"))
trait_dendro_pdf      <- file.path(out_dir, "Trait_Dendrogram_From_PC_Associations.pdf")
removed_non_numeric_file <- file.path(out_dir, "sample_info_non_numeric_columns_removed.txt")
removed_all_na_file      <- file.path(out_dir, "sample_info_all_na_numeric_columns_removed.txt")
removed_zero_var_file    <- file.path(out_dir, "sample_info_zero_variance_numeric_columns_removed.txt")

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
write_vector_file <- function(x, file) {
  x <- unique(as.character(x))
  if (length(x) == 0) x <- character(0)
  writeLines(x, con = file)
}

getPCs_generic <- function(meth,
                           mod = matrix(1, nrow = ncol(meth), ncol = 1),
                           save = TRUE,
                           file = "PCs.rds",
                           verbose = TRUE) {
  if (verbose) message("[getPCs_generic] Determining number of PCs using sva::num.sv")

  mod_noNA <- mod
  if (anyNA(mod_noNA)) {
    if (verbose) message("[getPCs_generic] NAs detected in model matrix; centering and replacing with 0")
    if (ncol(mod_noNA) > 1) {
      mod_noNA[, -1] <- scale(mod_noNA[, -1], center = TRUE, scale = TRUE)
    }
    mod_noNA[is.na(mod_noNA)] <- 0
  }

  n_pc <- sva::num.sv(meth, mod = mod_noNA, method = "be", seed = 5)
  if (verbose) message("[getPCs_generic] num.sv selected ", n_pc, " PCs as surrogate variables")

  meth_t <- t(meth)
  meth_t_centered <- scale(meth_t, center = TRUE, scale = FALSE)
  ss <- svd(meth_t_centered)

  PCs <- ss$u[, seq_len(n_pc), drop = FALSE]
  dimnames(PCs) <- list(rownames(meth_t), paste0("PC_", seq_len(n_pc)))

  var_explained     <- (ss$d^2) / sum(ss$d^2)
  cum_var_explained <- cumsum(var_explained)

  if (save) saveRDS(PCs, file = file)

  list(
    PCs               = PCs,
    singular_values   = ss$d,
    var_explained     = var_explained,
    cum_var_explained = cum_var_explained,
    n_pc              = n_pc
  )
}

label_sv_in_stats <- function(stats_df, n_sv) {
  stats_df$module <- as.character(stats_df$module)
  sv_pcs <- paste0("PC_", seq_len(n_sv))
  stats_df$module <- ifelse(stats_df$module %in% sv_pcs,
                            paste0(stats_df$module, " [SV]"),
                            stats_df$module)

  pc_nums <- suppressWarnings(as.integer(gsub("PC_([0-9]+).*", "\\1", unique(stats_df$module))))
  pc_nums <- sort(unique(pc_nums[!is.na(pc_nums)]))
  all_pcs <- paste0("PC_", pc_nums)
  new_levels <- ifelse(all_pcs %in% sv_pcs, paste0(all_pcs, " [SV]"), all_pcs)
  stats_df$module <- factor(stats_df$module, levels = new_levels)
  stats_df
}

# ------------------------------------------------------------
# Load methylation
# ------------------------------------------------------------
meth <- readRDS(opt$input_meth)
if (!(is.matrix(meth) || is.data.frame(meth))) stop("input_meth must contain a matrix-like object.")
meth <- as.matrix(meth)
if (!is.numeric(meth))       stop("Methylation matrix must be numeric.")
if (nrow(meth) < 1)          stop("Methylation matrix has zero rows.")
if (ncol(meth) < 2)          stop("Methylation matrix must have at least 2 samples.")
if (is.null(colnames(meth))) stop("Methylation matrix must have sample IDs in colnames.")
if (anyDuplicated(colnames(meth))) {
  dup_ids <- unique(colnames(meth)[duplicated(colnames(meth))])
  stop("Duplicate methylation sample IDs found. Example: ", paste(head(dup_ids, 10), collapse = ", "))
}
message("Loaded methylation matrix: ", nrow(meth), " regions x ", ncol(meth), " samples")

# ------------------------------------------------------------
# Load sample info
# ------------------------------------------------------------
sample_info <- readSampleInfo(
  file = opt$sample_info,
  sample_id_col = opt$sample_id_col,
  verbose = TRUE
)
if (nrow(sample_info) < 2) stop("Sample info must contain at least 2 samples after loading.")
if (is.null(rownames(sample_info))) stop("Sample info must have rownames after ID resolution.")
if (anyDuplicated(rownames(sample_info))) {
  dup_ids <- unique(rownames(sample_info)[duplicated(rownames(sample_info))])
  stop("Duplicate sample IDs found in sample info. Example: ", paste(head(dup_ids, 10), collapse = ", "))
}

# ------------------------------------------------------------
# Align samples and write alignment diagnostics
# ------------------------------------------------------------
meth_ids <- colnames(meth)
info_ids <- rownames(sample_info)
samples_only_in_meth <- setdiff(meth_ids, info_ids)
samples_only_in_info <- setdiff(info_ids, meth_ids)

alignment_dir <- file.path(out_dir, "sample_alignment")
dir.create(alignment_dir, recursive = TRUE, showWarnings = FALSE)
write_vector_file(samples_only_in_meth, file.path(alignment_dir, "samples_only_in_methylation.txt"))
write_vector_file(samples_only_in_info, file.path(alignment_dir, "samples_only_in_sample_info.txt"))

message("Samples only in methylation: ", length(samples_only_in_meth))
if (length(samples_only_in_meth) > 0) message("Examples only in methylation: ", paste(head(samples_only_in_meth, 10), collapse = ", "))
message("Samples only in sample_info: ", length(samples_only_in_info))
if (length(samples_only_in_info) > 0) message("Examples only in sample_info: ", paste(head(samples_only_in_info, 10), collapse = ", "))

common_samples <- intersect(meth_ids, info_ids)
if (length(common_samples) == 0) stop("No overlapping sample IDs between methylation matrix and sample info.")
if (length(common_samples) < 5) warning("Only ", length(common_samples), " overlapping samples found.")

meth <- meth[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]
message("After alignment: ", nrow(meth), " regions x ", ncol(meth), " samples")

# ------------------------------------------------------------
# Keep numeric traits only for PC-trait correlations
# ------------------------------------------------------------
is_numeric_atomic <- vapply(sample_info, function(x) is.numeric(x) && is.atomic(x), logical(1))
removed_non_numeric <- names(sample_info)[!is_numeric_atomic]
sample_info_num <- sample_info[, is_numeric_atomic, drop = FALSE]
write_vector_file(removed_non_numeric, removed_non_numeric_file)

if (ncol(sample_info_num) == 0) stop("No numeric trait columns available after filtering.")

all_na <- vapply(sample_info_num, function(x) all(is.na(x)), logical(1))
zero_var <- vapply(sample_info_num, function(x) {
  vals <- x[!is.na(x)]
  if (length(vals) <= 1) return(TRUE)
  length(unique(vals)) == 1
}, logical(1))

removed_all_na <- names(sample_info_num)[all_na]
removed_zero_var <- names(sample_info_num)[zero_var]
sample_info_num <- sample_info_num[, !(all_na | zero_var), drop = FALSE]

write_vector_file(removed_all_na, removed_all_na_file)
write_vector_file(removed_zero_var, removed_zero_var_file)

if (ncol(sample_info_num) == 0) stop("No usable numeric traits remain after removing all-NA and zero-variance columns.")
message("Numeric trait columns retained for diagnostics: ", ncol(sample_info_num))

# ------------------------------------------------------------
# Load and resolve protected / technical traits
# ------------------------------------------------------------
protected_traits_requested <- readTraitFile(opt$protected_traits_file, verbose = TRUE)
technical_traits_requested <- readTraitFile(opt$technical_traits_file, verbose = TRUE)

protected_resolved <- resolveTraits(
  requested_traits = protected_traits_requested,
  available_traits = colnames(sample_info_num),
  label = "protected_traits",
  verbose = TRUE
)

technical_resolved <- resolveTraits(
  requested_traits = technical_traits_requested,
  available_traits = colnames(sample_info_num),
  label = "technical_traits",
  verbose = TRUE
)

if (!is.null(opt$protected_traits_file)) {
  write_vector_file(protected_resolved$requested, file.path(out_dir, "traits_requested_protected.txt"))
  write_vector_file(protected_resolved$found,     file.path(out_dir, "traits_found_protected.txt"))
  write_vector_file(protected_resolved$missing,   file.path(out_dir, "traits_missing_protected.txt"))
}
if (!is.null(opt$technical_traits_file)) {
  write_vector_file(technical_resolved$requested, file.path(out_dir, "traits_requested_technical.txt"))
  write_vector_file(technical_resolved$found,     file.path(out_dir, "traits_found_technical.txt"))
  write_vector_file(technical_resolved$missing,   file.path(out_dir, "traits_missing_technical.txt"))
}

# ------------------------------------------------------------
# Determine selected trait mode for heatmaps
# ------------------------------------------------------------
selected_trait_mode <- opt$selected_trait_mode
if (is.null(selected_trait_mode)) {
  selected_trait_mode <- if (length(protected_resolved$found) > 0) "protected" else "none"
}
selected_trait_mode <- tolower(selected_trait_mode)

selected_traits <- character(0)
if (selected_trait_mode == "protected") {
  selected_traits <- protected_resolved$found
} else if (selected_trait_mode == "protected_and_technical") {
  selected_traits <- unique(c(protected_resolved$found, technical_resolved$found))
}

# ------------------------------------------------------------
# Build model matrix for num.sv
# ------------------------------------------------------------
old_na_action <- getOption("na.action")
options(na.action = "na.pass")
on.exit(options(na.action = old_na_action), add = TRUE)

if (opt$pc_model == "intercept_only") {
  mod_formula <- "~ 1"
  mod <- matrix(1, nrow = ncol(meth), ncol = 1)
  rownames(mod) <- colnames(meth)
  colnames(mod) <- "(Intercept)"
  message("PC model selected: intercept_only")
  message("Using intercept-only model for sva::num.sv")
} else if (opt$pc_model == "protected") {
  if (length(protected_resolved$found) == 0) {
    stop("--pc_model protected was requested, but no protected traits were found. ",
         "Either provide a valid --protected_traits_file or use --pc_model intercept_only.")
  }
  protected_df <- sample_info_num[, protected_resolved$found, drop = FALSE]
  mod_formula <- paste("~", paste(paste0("`", colnames(protected_df), "`"), collapse = " + ")) #mod_formula <- paste("~", paste(colnames(protected_df), collapse = " + "))
  mod <- model.matrix(as.formula(mod_formula), data = protected_df)
  message("PC model selected: protected")
  message("Using protected traits in SVA model matrix: ", paste(colnames(protected_df), collapse = ", "))
}

if (nrow(mod) != ncol(meth)) {
  stop("Model matrix row count (", nrow(mod), ") does not match methylation sample count (", ncol(meth), ").")
}

# ------------------------------------------------------------
# Derive PCs
# ------------------------------------------------------------
pc_out <- getPCs_generic(meth = meth, mod = mod, save = TRUE, file = pcs_rds, verbose = TRUE)
PCs <- pc_out$PCs
var_explained <- pc_out$var_explained
cum_var_explained <- pc_out$cum_var_explained
singular_values <- pc_out$singular_values
n_pc <- pc_out$n_pc

message("PC matrix dimensions: ", nrow(PCs), " samples x ", ncol(PCs), " PCs")
message("Surrogate variables (num.sv): PC_1 through PC_", n_pc)

# ------------------------------------------------------------
# Save PC metadata
# ------------------------------------------------------------
pc_meta <- data.frame(
  PC = paste0("PC_", seq_along(singular_values)),
  singular_value = singular_values,
  variance_explained = var_explained,
  cumulative_variance_explained = cum_var_explained,
  rank = seq_along(singular_values),
  is_surrogate_variable = seq_along(singular_values) <= n_pc,
  stringsAsFactors = FALSE
)
readr::write_tsv(pc_meta, pc_meta_tsv)

# ------------------------------------------------------------
# Variance explained plot
# ------------------------------------------------------------
n_plot <- min(opt$max_pcs_plot, length(var_explained))
var_df_plot <- data.frame(
  PC = seq_len(n_plot),
  variance_explained = var_explained[seq_len(n_plot)]
)
var_df_plot$is_sv <- var_df_plot$PC <= n_pc

p_var <- ggplot(var_df_plot, aes(x = PC, y = variance_explained)) +
  geom_line(color = "grey60", linewidth = 0.7) +
  geom_point(aes(color = is_sv), size = 2.5) +
  {if (n_pc > 0 && n_pc < nrow(var_df_plot))
    geom_vline(xintercept = n_pc + 0.5, linetype = "dashed", linewidth = 0.6, color = "firebrick")
  } +
  scale_color_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "grey40"),
    labels = c("TRUE" = paste0("Surrogate Variable (num.sv, n=", n_pc, ")"),
               "FALSE" = "PC (not selected as SV)"),
    name = NULL
  ) +
  scale_x_continuous(breaks = seq(1, nrow(var_df_plot), by = 2)) +
  theme_bw(base_size = opt$base_size) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(
    title = paste0("Variance Explained (Top ", opt$max_pcs_plot, " PCs | SVs selected by num.sv: ", n_pc, ")"),
    subtitle = paste0("PC model: ", pc_model_label, " | Red points = surrogate variables"),
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  )

ggsave(filename = var_plot_pdf, plot = p_var, width = opt$plot_width,
       height = opt$plot_height, dpi = opt$plot_dpi, units = "in", limitsize = FALSE)
message("Saved variance explained plot: ", var_plot_pdf)

# ------------------------------------------------------------
# PC-trait correlations
# ------------------------------------------------------------
PCtraitCor_pearson <- getMEtraitCor(
  PCs,
  colData = sample_info_num,
  corType = "pearson",
  file = pc_trait_pearson_tsv
)

PCtraitCor_bicor <- getMEtraitCor(
  PCs,
  colData = sample_info_num,
  corType = "bicor",
  file = pc_trait_bicor_tsv
)

primary_method <- opt$primary_cor_method
primary_stats <- if (primary_method == "bicor") PCtraitCor_bicor else PCtraitCor_pearson
primary_cor_column <- if (primary_method == "bicor") "bicor" else "cor"
primary_stats_sv <- label_sv_in_stats(primary_stats, n_sv = n_pc)

# ------------------------------------------------------------
# Selected traits heatmap
# ------------------------------------------------------------
if (length(selected_traits) > 0 && selected_trait_mode != "none") {
  plotPCTrait(
    MEtraitCor = primary_stats_sv,
    subsetTraits = selected_traits,
    file = heatmap_selected_pdf,
    width = opt$plot_width,
    height = opt$plot_height,
    dpi = opt$plot_dpi,
    base_size = opt$base_size,
    axis.text.size = opt$axis_text_size,
    verbose = TRUE
  )
} else {
  message("Skipping selected-traits heatmap: no selected traits available under selected_trait_mode = ", selected_trait_mode)
}

# ------------------------------------------------------------
# Top associations heatmap
# ------------------------------------------------------------
plotPCTrait(
  MEtraitCor = primary_stats_sv,
  topCells = opt$top_cell_n,
  file = heatmap_topcells_pdf,
  width = opt$plot_width,
  height = opt$plot_height,
  dpi = opt$plot_dpi,
  base_size = opt$base_size,
  axis.text.size = opt$axis_text_size,
  verbose = TRUE
)

# ------------------------------------------------------------
# Top traits overall heatmap
# ------------------------------------------------------------
top_traits_list <- names(sort(tapply(primary_stats_sv$p, primary_stats_sv$trait, min)))[
  seq_len(min(opt$top_trait_n, length(unique(primary_stats_sv$trait))))
]

plotPCTrait(
  MEtraitCor = primary_stats_sv[primary_stats_sv$trait %in% top_traits_list, ],
  file = heatmap_toptraits_pdf,
  width = opt$plot_width,
  height = opt$plot_height,
  dpi = opt$plot_dpi,
  base_size = opt$base_size,
  axis.text.size = opt$axis_text_size,
  verbose = TRUE
)

# ------------------------------------------------------------
# Trait dendrogram from PC associations
# ------------------------------------------------------------
plotTraitDendrogramFromPC(
  pc_trait_stats = primary_stats_sv,
  output_file = trait_dendro_pdf,
  cor_column = primary_cor_column,
  width = opt$plot_width,
  height = opt$plot_height,
  verbose = TRUE
)

# ------------------------------------------------------------
# Save run parameters
# ------------------------------------------------------------
run_params <- c(
  paste("project_root:", opt$project_root),
  paste("input_meth:", opt$input_meth),
  paste("sample_info:", opt$sample_info),
  paste("sample_id_col:", ifelse(is.null(opt$sample_id_col), "NULL", opt$sample_id_col)),
  paste("protected_traits_file:", ifelse(is.null(opt$protected_traits_file), "NULL", opt$protected_traits_file)),
  paste("technical_traits_file:", ifelse(is.null(opt$technical_traits_file), "NULL", opt$technical_traits_file)),
  paste("pc_model:", opt$pc_model),
  paste("pc_model_label:", pc_model_label),
  paste("threads:", opt$threads),
  paste("top_trait_n:", opt$top_trait_n),
  paste("top_cell_n:", opt$top_cell_n),
  paste("primary_cor_method:", primary_method),
  paste("selected_trait_mode:", selected_trait_mode),
  paste("max_pcs_plot:", opt$max_pcs_plot),
  paste("plot_width:", opt$plot_width),
  paste("plot_height:", opt$plot_height),
  paste("plot_dpi:", opt$plot_dpi),
  paste("base_size:", opt$base_size),
  paste("axis_text_size:", opt$axis_text_size),
  paste("cpg_label:", cpg_label),
  paste("region_label:", region_label),
  paste("model_formula:", mod_formula),
  paste("model_matrix_dimensions:", paste(dim(mod), collapse = " x ")),
  paste("n_samples_after_alignment:", ncol(meth)),
  paste("n_regions:", nrow(meth)),
  paste("n_numeric_traits_used:", ncol(sample_info_num)),
  paste("n_pcs_selected_as_svs_by_numsv:", n_pc),
  paste("variance_captured_by_selected_pcs_percent:", round(sum(var_explained[seq_len(n_pc)]) * 100, 2)),
  paste("n_protected_traits_requested:", length(protected_resolved$requested)),
  paste("n_protected_traits_found:", length(protected_resolved$found)),
  paste("n_technical_traits_requested:", length(technical_resolved$requested)),
  paste("n_technical_traits_found:", length(technical_resolved$found)),
  paste("samples_only_in_methylation:", length(samples_only_in_meth)),
  paste("samples_only_in_sample_info:", length(samples_only_in_info)),
  paste("date:", as.character(Sys.time()))
)

writeLines(run_params, con = file.path(out_dir, "run_parameters.txt"))

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------
message("Surrogate variables identified by num.sv: ", n_pc, " (PC_1 through PC_", n_pc, ")")
message("Consensus Script 04 complete: PC derivation and diagnostics finished")
message("PC model used for num.sv: ", opt$pc_model)
message("Output folder: ", out_dir)
message("------------------------------------------------------------")
message("[LATENT STRUCTURE SUMMARY]")
message("------------------------------------------------------------")
message("Samples: ", nrow(PCs))
message("Selected latent factors for adjustment (num.sv): ", n_pc)
message("")
message("Returned adjustment matrix: ", nrow(PCs), " x ", ncol(PCs), " (samples x selected PCs)")
message("")
message("Variance captured by selected PCs: ", round(sum(var_explained[seq_len(n_pc)]) * 100, 2), "%")
message("")
