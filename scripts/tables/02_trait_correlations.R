#!/usr/bin/env Rscript
# ============================================================
# 02_trait_correlations.R
# Descriptive statistics — Script 02: pairwise trait–trait correlations
#
# PURPOSE
#   Computes pairwise-complete correlations (Spearman by default) among the
#   traits listed in a trait file, corrects for multiple testing across the
#   unique trait pairs, writes a formatted Excel workbook, and draws annotated
#   heatmaps (one with unadjusted, one with adjusted significance stars).
#
# USAGE
#   pixi run Rscript 02_trait_correlations.R \
#     --input    data/metadata/Intergen173_05222026_publication_labels.xlsx \
#     --traits   config/presentation_sets/00_figure1_heatmap_primary_traits.txt \
#     --out_dir  results/descriptive_stats/02_trait_correlations \
#     --na_codes "F,N"
#
# REQUIRED INPUTS
#   --input     : trait table (.xlsx, .csv, or .tsv). Rows = samples,
#                 columns = traits, with one sample ID column. Column names
#                 are read verbatim (spaces, colons, slashes preserved).
#   --out_dir   : output directory; created if it does not exist.
#
# DATA SELECTION
#   --traits    : text file listing the traits to correlate, one per line.
#                 If omitted, every non-ID column is used. Format:
#                   - names must match --input column names exactly
#                   - lines starting with '#' are comments; blank lines ignored
#                   - section headers '# ── Group name ───' assign the traits
#                     below them to a group (heatmap divider lines, Group_1/
#                     Group_2 columns in pair tables, --pair_scope)
#                   - file order sets heatmap order
#                   - duplicate or unmatched names stop the run (unmatched
#                     names are reported with the 3 closest column names)
#   --sheet     : Excel sheet name or index; ignored for csv/tsv [default = 1]
#   --id_col    : sample ID column name; IDs must be unique and non-missing
#                 [default = first column]
#   --na_codes  : comma-separated cell codes to treat as missing in trait
#                 columns, e.g. "F" (failed assay in the CHDS dictionary).
#                 Exact, case-sensitive whole-cell match after trimming.
#                 Recoded counts are reported per trait. Any other
#                 non-numeric value stops the run [default = none]
#   --drop_constant : drop traits with < 2 distinct values (with a warning);
#                 if FALSE, stop instead [default = TRUE]
#
# STATISTICS
#   --method    : correlation method: spearman, pearson, or kendall
#                 [default = spearman]
#   --adjust    : multiple-testing correction; any p.adjust() method
#                 (holm, hochberg, hommel, bonferroni, BH, BY, fdr, none)
#                 [default = holm]
#   --pair_scope: which unique pairs form the adjustment family
#                   all            = every pair in the upper triangle
#                   between_groups = only pairs whose traits sit in different
#                                    trait-file sections; within-group pairs
#                                    get adjusted p = NA. Requires every trait
#                                    to be under a section header.
#                 [default = all]
#   --alpha     : threshold for the significant-pairs sheets and the p-value
#                 highlighting in the matrix sheets [default = 0.05]
#   --min_pairwise_n : pairs with fewer pairwise-complete samples than this
#                 get coefficient and p set to NA and are excluded from the
#                 adjustment family [default = 10]
#
# HEATMAP APPEARANCE
#   --star_thresholds : comma-separated p thresholds for 1, 2, 3, ... stars
#                 [default = 0.05,0.01,0.001,0.0001]
#   --show_group_lines: draw divider lines between trait-file sections
#                 [default = TRUE]
#   --width     : heatmap width in inches [default = auto, 0.30 × n_traits + 4]
#   --height    : heatmap height in inches [default = auto, 0.30 × n_traits + 3]
#   --axis_text_size : axis label font size in pt [default = 10]
#   --star_size : star text size in ggplot mm units [default = 3]
#
# OUTPUT NAMING
#   --prefix    : filename prefix for the workbook and heatmaps; use a
#                 distinct prefix per run variant (e.g. trait_correlations_
#                 between) to avoid overwriting [default = trait_correlations]
#
# OUTPUT
#   <out_dir>/
#       <prefix>.xlsx
#           Coefficients        : correlation matrix (colour scale)
#           P_unadjusted        : unadjusted p matrix (p < alpha highlighted)
#           P_<adjust>          : adjusted p matrix (p < alpha highlighted)
#           Pairwise_N          : pairwise-complete sample sizes
#           All_pairs           : one row per unique pair, sorted by p
#           Sig_unadj_p<alpha>  : pairs with unadjusted p < alpha
#           Sig_<adjust>_p<alpha>: pairs with adjusted p < alpha
#           Trait_summary       : N, missingness, recodes, median, range
#           Run_info            : parameters and counts for this run
#       heatmaps/
#           <prefix>_heatmap_unadjusted.pdf / .png
#           <prefix>_heatmap_<adjust>.pdf / .png
#       run_parameters.txt
#
# NOTES
#   - PDFs use cairo_pdf (vector text for Illustrator); PNGs are 600 DPI.
#   - helper.R must sit next to this script and define
#     setup_publication_font().
#   - Warnings from cor() ("standard deviation is zero", "NaNs produced")
#     come from pairs with too few overlapping samples; those pairs are
#     named in the console and set to NA via --min_pairwise_n.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
  library(psych)
  library(ggplot2)
})

# ── Locate and source helper.R ─────────────────────────────────────────────
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
}
script_dir  <- get_script_dir()
helper_path <- file.path(script_dir, "helper.R")
if (!file.exists(helper_path)) {
  stop("helper.R not found next to this script: ", helper_path)
}
source(helper_path)

# ── CLI arguments ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = paste("Path to the trait table (.xlsx, .csv, or .tsv).",
                           "Rows = samples, columns = traits. [required]")),
  make_option("--sheet", type = "character", default = "1",
              help = "Excel sheet name or index (ignored for csv/tsv). [default: %default]"),
  make_option("--id_col", type = "character", default = NULL,
              help = "Sample ID column name. [default: first column]"),
  make_option("--traits", type = "character", default = NULL,
              help = paste("Text file listing the traits to correlate (one per line,",
                           "'#' comments allowed). If omitted, all non-ID columns are used.")),
  make_option("--out_dir", type = "character", default = NULL,
              help = "Output directory (created if absent). [required]"),
  make_option("--prefix", type = "character", default = "trait_correlations",
              help = "Filename prefix for all outputs. [default: %default]"),
  make_option("--na_codes", type = "character", default = NULL,
              help = paste("Comma-separated cell codes to treat as missing in trait columns",
                           "(e.g. \"F\" for failed assay). Exact, case-sensitive match after",
                           "trimming whitespace. Recoded counts are reported per trait in",
                           "Trait_summary and run_parameters.txt. [default: none]")),
  make_option("--method", type = "character", default = "spearman",
              help = "Correlation method: spearman, pearson, or kendall. [default: %default]"),
  make_option("--adjust", type = "character", default = "holm",
              help = paste("Multiple-testing correction (any p.adjust method:",
                           paste(p.adjust.methods, collapse = ", "), "). [default: %default]")),
  make_option("--pair_scope", type = "character", default = "all",
              help = paste("Which unique pairs form the adjustment family:",
                           "'all' = every pair;",
                           "'between_groups' = only pairs whose traits are in different",
                           "trait-file sections (within-group pairs get adjusted p = NA).",
                           "[default: %default]")),
  make_option("--alpha", type = "double", default = 0.05,
              help = "Significance threshold for the significant-pairs sheets. [default: %default]"),
  make_option("--star_thresholds", type = "character", default = "0.05,0.01,0.001,0.0001",
              help = "Comma-separated p thresholds for 1..n heatmap stars. [default: %default]"),
  make_option("--min_pairwise_n", type = "integer", default = 10,
              help = paste("Minimum pairwise-complete N; pairs below this get",
                           "coefficient and p set to NA. [default: %default]")),
  make_option("--drop_constant", type = "logical", default = TRUE,
              help = paste("Drop traits with zero variance or fewer than 2 distinct",
                           "values (with a warning). If FALSE, stop instead. [default: %default]")),
  make_option("--show_group_lines", type = "logical", default = TRUE,
              help = "Draw divider lines between trait-file sections on heatmaps. [default: %default]"),
  make_option("--width", type = "double", default = NA,
              help = "Heatmap width in inches. [default: auto from number of traits]"),
  make_option("--height", type = "double", default = NA,
              help = "Heatmap height in inches. [default: auto from number of traits]"),
  make_option("--axis_text_size", type = "double", default = 10,
              help = "Axis label font size (pt). [default: %default]"),
  make_option("--star_size", type = "double", default = 3,
              help = "Star text size (ggplot mm units). [default: %default]")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "Pairwise trait–trait correlations with multiple-testing correction, Excel output, and heatmaps."
))

# ── Validation ─────────────────────────────────────────────────────────────
if (is.null(opt$input))   stop("--input is required.")
if (is.null(opt$out_dir)) stop("--out_dir is required.")
if (!file.exists(opt$input)) stop("Input file not found: ", opt$input)
if (!is.null(opt$traits) && !file.exists(opt$traits)) {
  stop("Trait file not found: ", opt$traits)
}
opt$method <- tolower(opt$method)
if (!opt$method %in% c("spearman", "pearson", "kendall")) {
  stop("--method must be one of: spearman, pearson, kendall.")
}
if (!opt$adjust %in% p.adjust.methods) {
  stop("--adjust must be one of: ", paste(p.adjust.methods, collapse = ", "))
}
if (!opt$pair_scope %in% c("all", "between_groups")) {
  stop("--pair_scope must be 'all' or 'between_groups'.")
}
if (!is.finite(opt$alpha) || opt$alpha <= 0 || opt$alpha >= 1) {
  stop("--alpha must be in (0, 1).")
}
star_thr <- suppressWarnings(as.numeric(strsplit(opt$star_thresholds, ",")[[1]]))
if (anyNA(star_thr) || any(star_thr <= 0 | star_thr >= 1)) {
  stop("--star_thresholds must be comma-separated numbers in (0, 1).")
}
star_thr <- sort(star_thr, decreasing = TRUE)
if (opt$min_pairwise_n < 3) stop("--min_pairwise_n must be >= 3.")

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
heatmap_dir <- file.path(opt$out_dir, "heatmaps")
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helpers ────────────────────────────────────────────────────────────────
na_strings <- c("NA", "N/A", "NaN", ".", "")

read_trait_table <- function(path, sheet) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "xlsx") {
    sheet_arg <- if (grepl("^[0-9]+$", sheet)) as.integer(sheet) else sheet
    # sep.names = " " keeps spaces in column names (openxlsx default replaces
    # them with '.'), so names like "F1 Black race" survive intact.
    openxlsx::read.xlsx(path, sheet = sheet_arg, check.names = FALSE,
                        sep.names = " ", na.strings = na_strings)
  } else if (ext == "csv") {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                    na.strings = na_strings)
  } else if (ext %in% c("tsv", "txt")) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                      na.strings = na_strings)
  } else {
    stop("Unsupported input extension '.", ext, "'. Use .xlsx, .csv, or .tsv.")
  }
}

parse_trait_file <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- trimws(sub("\r$", "", lines))
  # Section header: '#' followed by >= 2 box-drawing or hyphen characters
  header_re <- "^#\\s*[\u2500-]{2,}\\s*(.*?)\\s*[\u2500-]*$"
  current_group <- NA_character_
  traits <- character(0)
  groups <- character(0)
  for (ln in lines) {
    if (ln == "") next
    if (startsWith(ln, "#")) {
      if (grepl(header_re, ln, perl = TRUE)) {
        g <- trimws(sub(header_re, "\\1", ln, perl = TRUE))
        if (nzchar(g)) current_group <- g
      }
      next
    }
    ln <- trimws(sub("\\s+#.*$", "", ln))  # strip trailing inline comments
    if (ln == "") next
    traits <- c(traits, ln)
    groups <- c(groups, current_group)
  }
  if (length(traits) == 0) stop("No traits found in trait file: ", path)
  dups <- unique(traits[duplicated(traits)])
  if (length(dups) > 0) {
    stop("Duplicate traits in trait file: ", paste(shQuote(dups), collapse = ", "))
  }
  list(traits = traits, groups = groups)
}

suggest_names <- function(missing, available) {
  vapply(missing, function(m) {
    d <- utils::adist(tolower(m), tolower(available))[1, ]
    best <- available[order(d)][1:min(3, length(available))]
    paste0("  - ", shQuote(m), "  (closest: ", paste(shQuote(best), collapse = ", "), ")")
  }, character(1))
}

coerce_numeric <- function(df, codes = character(0)) {
  bad <- list()
  recoded <- setNames(integer(ncol(df)), names(df))
  for (nm in names(df)) {
    x <- df[[nm]]
    if (is.factor(x))  x <- as.character(x)
    if (is.logical(x)) x <- as.integer(x)
    if (length(codes) > 0 && is.character(x)) {
      hit <- !is.na(x) & trimws(x) %in% codes
      recoded[nm] <- sum(hit)
      x[hit] <- NA
    }
    xn <- suppressWarnings(as.numeric(x))
    lost <- !is.na(x) & is.na(xn)
    if (any(lost)) bad[[nm]] <- utils::head(unique(as.character(x[lost])), 5)
    df[[nm]] <- xn
  }
  list(df = df, bad = bad, recoded = recoded)
}

p_to_stars <- function(p, thr) {
  out <- rep("", length(p))
  for (i in seq_along(thr)) out[!is.na(p) & p < thr[i]] <- strrep("*", i)
  out
}

# ── Load data ──────────────────────────────────────────────────────────────
message("Reading input: ", opt$input)
raw <- read_trait_table(opt$input, opt$sheet)

id_col <- if (is.null(opt$id_col)) names(raw)[1] else opt$id_col
if (!id_col %in% names(raw)) stop("ID column not found: ", shQuote(id_col))
ids <- as.character(raw[[id_col]])
if (anyNA(ids) || any(ids == "")) stop("ID column ", shQuote(id_col), " has missing values.")
if (anyDuplicated(ids)) {
  stop("Duplicate sample IDs in ", shQuote(id_col), ": ",
       paste(utils::head(unique(ids[duplicated(ids)]), 10), collapse = ", "))
}
available <- setdiff(names(raw), id_col)

# ── Resolve trait list ─────────────────────────────────────────────────────
if (is.null(opt$traits)) {
  traits <- available
  groups <- rep(NA_character_, length(traits))
  message("No --traits file given; using all ", length(traits), " non-ID columns.")
} else {
  tf <- parse_trait_file(opt$traits)
  traits <- tf$traits
  groups <- tf$groups
  missing <- setdiff(traits, available)
  if (length(missing) > 0) {
    stop("Traits in ", basename(opt$traits), " not found in input columns:\n",
         paste(suggest_names(missing, available), collapse = "\n"))
  }
  message("Using ", length(traits), " traits from ", basename(opt$traits),
          if (all(!is.na(groups))) paste0(" in ", length(unique(groups)), " groups") else "")
}

has_groups <- all(!is.na(groups)) && length(unique(groups)) > 1
if (opt$pair_scope == "between_groups" && !has_groups) {
  stop("--pair_scope between_groups requires every trait to sit under a ",
       "'# ── Section ──' header in the trait file, with at least two sections.")
}

# ── Coerce to numeric ──────────────────────────────────────────────────────
na_codes <- if (is.null(opt$na_codes)) character(0) else
  unique(trimws(strsplit(opt$na_codes, ",")[[1]]))
na_codes <- na_codes[nzchar(na_codes)]

co <- coerce_numeric(raw[, traits, drop = FALSE], na_codes)
if (length(co$bad) > 0) {
  msg <- vapply(names(co$bad), function(nm) {
    paste0("  - ", shQuote(nm), ": e.g. ", paste(shQuote(co$bad[[nm]]), collapse = ", "))
  }, character(1))
  stop("Non-numeric values found. Either recode them upstream, or, if they are ",
       "documented missing-value codes, pass them via --na_codes:\n",
       paste(msg, collapse = "\n"))
}
recoded <- co$recoded[co$recoded > 0]
if (length(recoded) > 0) {
  message("Recoded to NA via --na_codes (", paste(na_codes, collapse = ","), "):\n",
          paste0("  - ", names(recoded), ": ", recoded, collapse = "\n"))
}
dat <- co$df
rownames(dat) <- ids

# ── Drop uninformative traits ──────────────────────────────────────────────
n_distinct <- vapply(dat, function(x) length(unique(x[!is.na(x)])), integer(1))
constant <- names(n_distinct)[n_distinct < 2]
if (length(constant) > 0) {
  if (!opt$drop_constant) {
    stop("Traits with < 2 distinct values: ", paste(shQuote(constant), collapse = ", "))
  }
  warning("Dropping traits with < 2 distinct values: ",
          paste(shQuote(constant), collapse = ", "), call. = FALSE)
  keep   <- !traits %in% constant
  traits <- traits[keep]
  groups <- groups[keep]
  dat    <- dat[, traits, drop = FALSE]
}
if (length(traits) < 2) stop("Fewer than 2 usable traits remain.")
k <- length(traits)

# ── Correlations ───────────────────────────────────────────────────────────
message("Computing ", opt$method, " correlations for ", k, " traits (",
        k * (k - 1) / 2, " unique pairs)...")
mat <- as.matrix(dat)
ct  <- psych::corr.test(mat, use = "pairwise", method = opt$method,
                        adjust = "none", ci = FALSE)
r <- ct$r
p_unadj <- ct$p
p_unadj[upper.tri(p_unadj)] <- t(p_unadj)[upper.tri(p_unadj)]  # symmetric
diag(p_unadj) <- NA

# Exact pairwise-complete N (corr.test returns a scalar when nothing is missing)
n_mat <- crossprod(!is.na(mat))
dimnames(n_mat) <- dimnames(r)

low_n <- n_mat < opt$min_pairwise_n
diag(low_n) <- FALSE
if (any(low_n)) {
  ln_idx <- which(low_n & upper.tri(low_n), arr.ind = TRUE)
  message("Setting ", nrow(ln_idx), " pair(s) with N < ", opt$min_pairwise_n,
          " to NA (zero-SD / NaN warnings from cor() above usually come from these):\n",
          paste0("  - ", traits[ln_idx[, 1]], " x ", traits[ln_idx[, 2]],
                 " (N = ", n_mat[ln_idx], ")", collapse = "\n"))
  r[low_n] <- NA
  p_unadj[low_n] <- NA
}

# ── Multiple-testing correction across unique pairs ────────────────────────
in_family <- upper.tri(r)
if (opt$pair_scope == "between_groups") {
  in_family <- in_family & outer(groups, groups, "!=")
}
fam_idx <- which(in_family & !is.na(p_unadj))
p_adj <- matrix(NA_real_, k, k, dimnames = dimnames(r))
p_adj[fam_idx] <- p.adjust(p_unadj[fam_idx], method = opt$adjust)
p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]
n_tests <- length(fam_idx)
message("Adjustment family (", opt$adjust, ", scope = ", opt$pair_scope, "): ",
        n_tests, " tests.")

# ── Long-format pair tables ────────────────────────────────────────────────
adj_col <- paste0("P_", opt$adjust)
ut <- which(upper.tri(r), arr.ind = TRUE)
pairs <- data.frame(
  Trait_1      = traits[ut[, 1]],
  Group_1      = groups[ut[, 1]],
  Trait_2      = traits[ut[, 2]],
  Group_2      = groups[ut[, 2]],
  N            = n_mat[ut],
  Coefficient  = r[ut],
  P_unadjusted = p_unadj[ut],
  P_adjusted   = p_adj[ut],
  In_adjustment_family = in_family[ut],
  check.names = FALSE, stringsAsFactors = FALSE
)
names(pairs)[names(pairs) == "P_adjusted"] <- adj_col
if (!has_groups) pairs$Group_1 <- pairs$Group_2 <- NULL
pairs <- pairs[order(pairs$P_unadjusted, na.last = TRUE), ]
rownames(pairs) <- NULL

sig_unadj <- pairs[!is.na(pairs$P_unadjusted) & pairs$P_unadjusted < opt$alpha, ]
sig_adj   <- pairs[!is.na(pairs[[adj_col]])   & pairs[[adj_col]]   < opt$alpha, ]

# ── Per-trait summary ──────────────────────────────────────────────────────
trait_summary <- data.frame(
  Trait       = traits,
  Group       = groups,
  N_nonmissing = vapply(dat, function(x) sum(!is.na(x)), integer(1)),
  N_missing    = vapply(dat, function(x) sum(is.na(x)), integer(1)),
  N_distinct   = vapply(dat, function(x) length(unique(x[!is.na(x)])), integer(1)),
  Median       = vapply(dat, function(x) stats::median(x, na.rm = TRUE), numeric(1)),
  Min          = vapply(dat, function(x) min(x, na.rm = TRUE), numeric(1)),
  Max          = vapply(dat, function(x) max(x, na.rm = TRUE), numeric(1)),
  check.names = FALSE, stringsAsFactors = FALSE
)
trait_summary$Pct_missing <- round(100 * trait_summary$N_missing / nrow(dat), 1)
trait_summary$N_recoded_na_codes <- unname(co$recoded[traits])
if (!has_groups) trait_summary$Group <- NULL
rownames(trait_summary) <- NULL

# ── Excel workbook ─────────────────────────────────────────────────────────
xlsx_path <- file.path(opt$out_dir, paste0(opt$prefix, ".xlsx"))
message("Writing workbook: ", xlsx_path)

wb       <- createWorkbook()
hdr      <- createStyle(textDecoration = "bold")
pos_sty  <- createStyle(fontColour = "#B2182B")
neg_sty  <- createStyle(fontColour = "#2166AC")
sig_sty  <- createStyle(textDecoration = "bold", bgFill = "#FFF2CC")
num3     <- createStyle(numFmt = "0.000")
sci      <- createStyle(numFmt = "0.00E+00")

add_matrix_sheet <- function(name, m, style = NULL) {
  addWorksheet(wb, name)
  # keepNA writes "NA" text rather than blanks: Excel treats blank cells as 0
  # in conditional rules, which would falsely flag the diagonal as p < alpha.
  writeData(wb, name, data.frame(m, check.names = FALSE), rowNames = TRUE,
            headerStyle = hdr, keepNA = TRUE, na.string = "NA")
  freezePane(wb, name, firstRow = TRUE, firstCol = TRUE)
  setColWidths(wb, name, cols = 1, widths = "auto")
  if (!is.null(style)) {
    addStyle(wb, name, style, rows = 2:(k + 1), cols = 2:(k + 1),
             gridExpand = TRUE, stack = TRUE)
  }
}

add_table_sheet <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, rowNames = FALSE, headerStyle = hdr)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
  if (nrow(df) == 0) return(invisible())
  rows <- 2:(nrow(df) + 1)
  # Locate columns by name so formatting never drifts to the wrong column
  cc <- match("Coefficient", names(df))
  if (!is.na(cc)) {
    addStyle(wb, name, num3, rows = rows, cols = cc, stack = TRUE)
    conditionalFormatting(wb, name, cols = cc, rows = rows, rule = "<0", style = neg_sty)
    conditionalFormatting(wb, name, cols = cc, rows = rows, rule = ">0", style = pos_sty)
  }
  for (pc in c("P_unadjusted", adj_col)) {
    j <- match(pc, names(df))
    if (!is.na(j)) addStyle(wb, name, sci, rows = rows, cols = j, stack = TRUE)
  }
}

# Matrices
add_matrix_sheet("Coefficients", r, num3)
conditionalFormatting(wb, "Coefficients", cols = 2:(k + 1), rows = 2:(k + 1),
                      style = c("#2166AC", "#FFFFFF", "#B2182B"),
                      rule = c(-1, 0, 1), type = "colourScale")
add_matrix_sheet("P_unadjusted", p_unadj, sci)
conditionalFormatting(wb, "P_unadjusted", cols = 2:(k + 1), rows = 2:(k + 1),
                      rule = paste0("<", opt$alpha), style = sig_sty)
add_matrix_sheet(adj_col, p_adj, sci)
conditionalFormatting(wb, adj_col, cols = 2:(k + 1), rows = 2:(k + 1),
                      rule = paste0("<", opt$alpha), style = sig_sty)
add_matrix_sheet("Pairwise_N", n_mat)

# Long tables
add_table_sheet("All_pairs", pairs)
add_table_sheet(paste0("Sig_unadj_p", opt$alpha), sig_unadj)
add_table_sheet(paste0("Sig_", opt$adjust, "_p", opt$alpha), sig_adj)
add_table_sheet("Trait_summary", trait_summary)

# Run info
run_info <- data.frame(
  Parameter = c("Input", "Sheet", "ID column", "Trait file", "Method",
                "Adjustment", "Pair scope", "Tests in adjustment family",
                "Alpha", "Min pairwise N", "NA codes", "Cells recoded to NA",
                "Samples", "Traits used",
                "Traits dropped (constant)", "Run at"),
  Value = c(normalizePath(opt$input), opt$sheet, id_col,
            if (is.null(opt$traits)) "(all columns)" else normalizePath(opt$traits),
            opt$method, opt$adjust, opt$pair_scope, n_tests, opt$alpha,
            opt$min_pairwise_n,
            if (length(na_codes)) paste(na_codes, collapse = ", ") else "none",
            if (length(recoded)) paste0(names(recoded), " (", recoded, ")", collapse = "; ") else "none",
            nrow(dat), k,
            if (length(constant)) paste(constant, collapse = "; ") else "none",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  stringsAsFactors = FALSE
)
add_table_sheet("Run_info", run_info)

saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# ── Heatmaps ───────────────────────────────────────────────────────────────
if (!exists("setup_publication_font", mode = "function")) {
  stop(
    "helper.R must define setup_publication_font(). ",
    "Add the reusable publication-font helper before running this script."
  )
}
font_family <- setup_publication_font(
  candidates = c("Arial", "Helvetica"),
  strict = TRUE
)

fill_scale <- if (requireNamespace("WGCNA", quietly = TRUE)) {
  scale_fill_gradientn(colors = WGCNA::blueWhiteRed(100, gamma = 0.9),
                       limits = c(-1, 1), na.value = "grey85",
                       name = paste0(tools::toTitleCase(opt$method), "\ncorrelation"))
} else {
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-1, 1), na.value = "grey85",
                       name = paste0(tools::toTitleCase(opt$method), "\ncorrelation"))
}

plot_heatmap <- function(r, p, star_label) {
  df <- data.frame(
    x     = factor(rep(traits, times = k), levels = traits),
    y     = factor(rep(traits, each  = k), levels = rev(traits)),
    rho   = as.vector(r),
    stars = p_to_stars(as.vector(p), star_thr)
  )
  thr_txt <- paste(vapply(seq_along(star_thr), function(i) {
    paste0(strrep("*", i), " p < ", format(star_thr[i], scientific = FALSE))
  }, character(1)), collapse = ";  ")

  g <- ggplot(df, aes(x, y, fill = rho)) +
    geom_tile() +
    geom_text(aes(label = stars), size = opt$star_size, family = font_family,
              color = "black", vjust = 0.75) +
    fill_scale +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, caption = paste0(thr_txt, "  (", star_label, ")")) +
    theme(
      text             = element_text(color = "black", family = font_family,
                                      size = opt$axis_text_size),
      axis.text        = element_text(color = "black", size = opt$axis_text_size),
      axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.ticks       = element_line(color = "black"),
      panel.background = element_blank(),
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.justification = "top",
      plot.caption     = element_text(hjust = 0, size = opt$axis_text_size * 0.8),
      plot.margin      = margin(5, 5, 5, 5)
    )

  if (opt$show_group_lines && has_groups) {
    b <- which(utils::head(groups, -1) != utils::tail(groups, -1))
    g <- g +
      geom_vline(xintercept = b + 0.5, linewidth = 0.6, color = "black") +
      geom_hline(yintercept = k - b + 0.5, linewidth = 0.6, color = "black")
  }
  g
}

save_plot <- function(g, base, width, height) {
  ggsave(paste0(base, ".pdf"), g, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  has_showtext <- requireNamespace("showtext", quietly = TRUE)
  if (has_showtext) showtext::showtext_opts(dpi = 600)
  ggsave(paste0(base, ".png"), g, width = width, height = height,
         dpi = 600, limitsize = FALSE)
  if (has_showtext) showtext::showtext_opts(dpi = 96)
  message("Saved: ", base, ".pdf / .png")
}

w <- if (is.na(opt$width))  max(7, 0.30 * k + 4) else opt$width
h <- if (is.na(opt$height)) max(7, 0.30 * k + 3) else opt$height

save_plot(plot_heatmap(r, p_unadj, "unadjusted"),
          file.path(heatmap_dir, paste0(opt$prefix, "_heatmap_unadjusted")), w, h)
adj_label <- paste0(opt$adjust, "-adjusted",
                    if (opt$pair_scope == "between_groups") ", between-group pairs only" else "")
save_plot(plot_heatmap(r, p_adj, adj_label),
          file.path(heatmap_dir, paste0(opt$prefix, "_heatmap_", opt$adjust)), w, h)

# ── run_parameters.txt ─────────────────────────────────────────────────────
opt_log <- opt[setdiff(names(opt), "help")]
log_lines <- c(
  "02_trait_correlations.R — run parameters",
  paste("Run at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Script dir:", script_dir),
  "",
  "## Arguments",
  sprintf("%s = %s", names(opt_log),
          vapply(opt_log, function(v) paste(v, collapse = ","), character(1))),
  "",
  "## Data",
  paste("Samples:", nrow(dat)),
  paste("NA codes:", if (length(na_codes)) paste(na_codes, collapse = ", ") else "none"),
  paste("Cells recoded to NA:",
        if (length(recoded)) paste0(names(recoded), " (", recoded, ")", collapse = "; ") else "none"),
  paste("Traits requested:", k + length(constant)),
  paste("Traits used:", k),
  paste("Traits dropped (constant):",
        if (length(constant)) paste(constant, collapse = "; ") else "none"),
  paste("Unique pairs:", k * (k - 1) / 2),
  paste("Pairs set to NA (N <", opt$min_pairwise_n, "):", sum(low_n[upper.tri(low_n)])),
  paste("Tests in adjustment family:", n_tests),
  paste("Significant pairs, unadjusted p <", opt$alpha, ":", nrow(sig_unadj)),
  paste("Significant pairs,", opt$adjust, "p <", opt$alpha, ":", nrow(sig_adj)),
  "",
  "## Software",
  R.version.string,
  paste0("psych ", utils::packageVersion("psych"),
         "; openxlsx ", utils::packageVersion("openxlsx"),
         "; ggplot2 ", utils::packageVersion("ggplot2"))
)
writeLines(log_lines, file.path(opt$out_dir, "run_parameters.txt"))

message("Done. ", nrow(sig_unadj), " pairs significant unadjusted; ",
        nrow(sig_adj), " after ", opt$adjust, " (family = ", n_tests, " tests).")