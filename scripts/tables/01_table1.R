#!/usr/bin/env Rscript
# ============================================================
# 01_table1.R
# Descriptive statistics — Script 01: Table 1 sample characteristics
#
# PURPOSE
#   Builds Table 1 from the publication-label metadata, using a trait file
#   for the rows so every label is spelled exactly as in the figures.
#   Rows are summarised overall and within each stratum (race, obesity,
#   menarche timing by default), with a between-group test per stratum.
#   Writes a formatted Excel workbook and a vector PDF.
#
# USAGE
#   pixi run Rscript scripts/tables/01_table1.R \
#     --input    data/metadata/Intergen173_05222026_publication_labels.xlsx \
#     --traits   config/presentation_sets/table1_traits.txt \
#     --out_dir  results/descriptive_stats/01_table1 \
#     --na_codes "F,N" \
#     --mean_sd_traits "F1 age (years),F1 age at menarche"
#
# REQUIRED INPUTS
#   --input     : trait table (.xlsx, .csv, or .tsv), e.g. the output of
#                 rename_metadata_publication_labels.R. Rows = samples,
#                 one sample ID column. Column names are read verbatim.
#   --traits    : text file listing the Table 1 rows, one per line.
#                   - names must match --input column names exactly;
#                     unmatched names stop the run and print the 3
#                     closest column names
#                   - '#' lines are comments; blank lines ignored
#                   - section headers '# ── Group name ───' become the
#                     row-group headers in the table
#                   - file order sets row order
#   --out_dir   : output directory; created if it does not exist.
#
# DATA SELECTION
#   --sheet     : Excel sheet name or index; ignored for csv/tsv [default = 1]
#   --id_col    : sample ID column; IDs must be unique and non-missing
#                 [default = first column]
#   --na_codes  : comma-separated cell codes to treat as missing, e.g. "F"
#                 (failed assay in the CHDS dictionary). Exact, case-
#                 sensitive whole-cell match after trimming. Recoded counts
#                 are reported in Variable_info. Any other non-numeric
#                 value stops the run [default = none]
#
# STRATIFICATION
#   --strata    : comma-separated strata to include, in column order. Any of
#                 race, obesity, menarche [default = race,obesity,menarche]
#   --race_col  : 0/1 column for race; 1 = Black [default = F1 Black race]
#   --obesity_col : 0/1 column for obesity [default = F1 obesity]
#   --early_menarche_col / --late_menarche_col : 0/1 indicator columns used
#                 to derive menarche timing (early = 1 → Early; late = 1 →
#                 Late; both 0 → Typical; otherwise missing). A participant
#                 flagged both early and late stops the run.
#                 [defaults = F1 early menarche / F1 late menarche]
#   --menarche_age_col : continuous age at menarche; treated as a defining
#                 variable of the menarche stratum [default = F1 age at menarche]
#   --menarche_labels : comma-separated display labels for Early, Typical,
#                 Late [default = "Early (<12 y),Typical (12–14 y),Late (>14 y)"]
#
#   A stratum's defining variables (e.g. "F1 Black race" in the race stratum)
#   are shown as "—" in that stratum's columns: their split is true by
#   construction and the test would be trivially significant.
#
# STATISTICS
#   Variable type is detected from the data: columns whose non-missing
#   values are all 0/1 are binary, everything else is continuous.
#   --cont_stat : continuous summary: median_iqr (median [Q1–Q3]) or
#                 mean_sd (mean (SD)) [default = median_iqr]
#   --mean_sd_traits : comma-separated traits summarised as mean (SD)
#                 regardless of --cont_stat, e.g. age [default = none]
#   --binary_test : fisher (exact) or chisq (Pearson, no continuity
#                 correction) [default = fisher]. Continuous variables use
#                 Wilcoxon rank-sum (2 groups) or Kruskal–Wallis (3+).
#   --alpha     : threshold for bolding p-values when --highlight_sig is
#                 TRUE [default = 0.05]
#
# FORMATTING
#   --digits    : significant digits for continuous summaries [default = 3]
#   --highlight_sig : bold p < alpha in Excel and PDF [default = TRUE]
#   --title     : table title [default = "Table 1. Characteristics of the
#                 study sample"]
#   --prefix    : filename prefix for outputs [default = table1]
#   --verify    : print a data verification report to the console before
#                 building the table (types, missingness, strata counts)
#                 [default = FALSE]
#
# OUTPUT
#   <out_dir>/
#       <prefix>.xlsx
#           Table1          : combined table (Overall + every stratum)
#           <stratum>       : one sheet per stratum
#           Variable_info   : type, summary, N, missing, recoded per trait
#           Run_info        : parameters and counts for this run
#       <prefix>.pdf        : page 1 combined table, then one page per stratum
#       run_parameters.txt
#
# NOTES
#   - helper.R must sit next to this script and define read_trait_table(),
#     parse_trait_file(), suggest_names(), coerce_numeric() and
#     setup_publication_font().
#   - PDF uses cairo_pdf (vector text for Illustrator).
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
  library(grid)
  library(gridExtra)
  library(gtable)
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

needed_helpers <- c("read_trait_table", "parse_trait_file", "suggest_names",
                    "coerce_numeric", "setup_publication_font")
missing_helpers <- needed_helpers[!vapply(needed_helpers, exists, logical(1),
                                          mode = "function")]
if (length(missing_helpers) > 0) {
  stop("helper.R is missing: ", paste(missing_helpers, collapse = ", "),
       ". Append helper_additions_trait_io.R to helper.R.")
}

# ── CLI arguments ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Path to the trait table (.xlsx, .csv, or .tsv). [required]"),
  make_option("--traits", type = "character", default = NULL,
              help = "Text file listing Table 1 rows, with '# ── Section ──' headers. [required]"),
  make_option("--out_dir", type = "character", default = NULL,
              help = "Output directory (created if absent). [required]"),
  make_option("--sheet", type = "character", default = "1",
              help = "Excel sheet name or index (ignored for csv/tsv). [default: %default]"),
  make_option("--id_col", type = "character", default = NULL,
              help = "Sample ID column name. [default: first column]"),
  make_option("--na_codes", type = "character", default = NULL,
              help = paste("Comma-separated cell codes to treat as missing (e.g. \"F\" for",
                           "failed assay). Exact, case-sensitive match after trimming. [default: none]")),
  make_option("--strata", type = "character", default = "race,obesity,menarche",
              help = "Comma-separated strata, in column order: race, obesity, menarche. [default: %default]"),
  make_option("--race_col", type = "character", default = "F1 Black race",
              help = "0/1 race column; 1 = Black. [default: %default]"),
  make_option("--obesity_col", type = "character", default = "F1 obesity",
              help = "0/1 obesity column; 1 = obese. [default: %default]"),
  make_option("--early_menarche_col", type = "character", default = "F1 early menarche",
              help = "0/1 early-menarche indicator. [default: %default]"),
  make_option("--late_menarche_col", type = "character", default = "F1 late menarche",
              help = "0/1 late-menarche indicator. [default: %default]"),
  make_option("--menarche_age_col", type = "character", default = "F1 age at menarche",
              help = "Continuous age at menarche (defining variable of the menarche stratum). [default: %default]"),
  make_option("--menarche_labels", type = "character",
              default = "Early (<12 y),Typical (12–14 y),Late (>14 y)",
              help = "Display labels for Early, Typical, Late. [default: %default]"),
  make_option("--cont_stat", type = "character", default = "median_iqr",
              help = "Continuous summary: median_iqr or mean_sd. [default: %default]"),
  make_option("--mean_sd_traits", type = "character", default = NULL,
              help = "Comma-separated traits shown as mean (SD) regardless of --cont_stat. [default: none]"),
  make_option("--binary_test", type = "character", default = "fisher",
              help = "Test for binary variables: fisher or chisq. [default: %default]"),
  make_option("--alpha", type = "double", default = 0.05,
              help = "Threshold for bolding p-values. [default: %default]"),
  make_option("--digits", type = "integer", default = 3,
              help = "Significant digits for continuous summaries. [default: %default]"),
  make_option("--highlight_sig", type = "logical", default = TRUE,
              help = "Bold p < alpha in Excel and PDF. [default: %default]"),
  make_option("--title", type = "character",
              default = "Table 1. Characteristics of the study sample",
              help = "Table title. [default: %default]"),
  make_option("--prefix", type = "character", default = "table1",
              help = "Filename prefix for all outputs. [default: %default]"),
  make_option("--verify", action = "store_true", default = FALSE,
              help = "Print a data verification report before building the table. [default: %default]")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "Table 1 sample characteristics from a trait file, stratified, with Excel and PDF output."
))

# ── Validation ─────────────────────────────────────────────────────────────
if (is.null(opt$input))   stop("--input is required.")
if (is.null(opt$traits))  stop("--traits is required.")
if (is.null(opt$out_dir)) stop("--out_dir is required.")
if (!file.exists(opt$input))  stop("Input file not found: ", opt$input)
if (!file.exists(opt$traits)) stop("Trait file not found: ", opt$traits)
if (!opt$cont_stat %in% c("median_iqr", "mean_sd")) {
  stop("--cont_stat must be 'median_iqr' or 'mean_sd'.")
}
if (!opt$binary_test %in% c("fisher", "chisq")) {
  stop("--binary_test must be 'fisher' or 'chisq'.")
}
if (!is.finite(opt$alpha) || opt$alpha <= 0 || opt$alpha >= 1) {
  stop("--alpha must be in (0, 1).")
}
if (opt$digits < 1) stop("--digits must be >= 1.")

split_csv <- function(s) {
  if (is.null(s)) return(character(0))
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}

strata_keys <- tolower(split_csv(opt$strata))
bad_strata  <- setdiff(strata_keys, c("race", "obesity", "menarche"))
if (length(bad_strata) > 0) {
  stop("Unknown --strata: ", paste(bad_strata, collapse = ", "),
       ". Use any of race, obesity, menarche.")
}
if (anyDuplicated(strata_keys)) stop("--strata lists a stratum twice.")

menarche_labels <- split_csv(opt$menarche_labels)
if (length(menarche_labels) != 3) {
  stop("--menarche_labels needs exactly 3 labels (Early, Typical, Late).")
}

na_codes       <- unique(split_csv(opt$na_codes))
mean_sd_traits <- split_csv(opt$mean_sd_traits)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

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
n_total   <- nrow(raw)

# ── Resolve trait list ─────────────────────────────────────────────────────
tf     <- parse_trait_file(opt$traits)
traits <- tf$traits
groups <- tf$groups

missing <- setdiff(traits, available)
if (length(missing) > 0) {
  stop("Traits in ", basename(opt$traits), " not found in input columns:\n",
       paste(suggest_names(missing, available), collapse = "\n"))
}
unknown_msd <- setdiff(mean_sd_traits, traits)
if (length(unknown_msd) > 0) {
  stop("--mean_sd_traits not in the trait file: ",
       paste(shQuote(unknown_msd), collapse = ", "))
}

# Strata columns must exist even if they are not Table 1 rows
strata_cols <- c(
  if ("race" %in% strata_keys)     opt$race_col,
  if ("obesity" %in% strata_keys)  opt$obesity_col,
  if ("menarche" %in% strata_keys) c(opt$early_menarche_col, opt$late_menarche_col)
)
missing_strata <- setdiff(strata_cols, available)
if (length(missing_strata) > 0) {
  stop("Stratum columns not found in input:\n",
       paste(suggest_names(missing_strata, available), collapse = "\n"))
}

# ── Coerce to numeric ──────────────────────────────────────────────────────
use_cols <- unique(c(traits, strata_cols))
co <- coerce_numeric(raw[, use_cols, drop = FALSE], na_codes)
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

# ── Variable types ─────────────────────────────────────────────────────────
is_binary01 <- function(x) {
  v <- unique(x[!is.na(x)])
  length(v) > 0 && all(v %in% c(0, 1))
}
all_missing <- traits[vapply(traits, function(t) all(is.na(dat[[t]])), logical(1))]
if (length(all_missing) > 0) {
  stop("Traits with no non-missing values: ",
       paste(shQuote(all_missing), collapse = ", "))
}
trait_type <- stats::setNames(
  ifelse(vapply(traits, function(t) is_binary01(dat[[t]]), logical(1)),
         "binary", "continuous"),
  traits)
bin_msd <- mean_sd_traits[trait_type[mean_sd_traits] == "binary"]
if (length(bin_msd) > 0) {
  stop("--mean_sd_traits includes binary variables: ",
       paste(shQuote(bin_msd), collapse = ", "))
}

# ── Strata ─────────────────────────────────────────────────────────────────
check_01 <- function(col) {
  x   <- dat[[col]]
  bad <- setdiff(unique(x[!is.na(x)]), c(0, 1))
  if (length(bad) > 0) {
    stop("Stratum column ", shQuote(col), " must be coded 0/1; found: ",
         paste(utils::head(bad, 5), collapse = ", "))
  }
  x
}

make_binary_stratum <- function(key, label, sheet, col, lab0, lab1) {
  x <- check_01(col)
  g <- factor(ifelse(is.na(x), NA, ifelse(x == 1, lab1, lab0)),
              levels = c(lab0, lab1))
  list(key = key, label = label, sheet = sheet, group = g, defining = col,
       definition = sprintf("%s: %s (%s = 1) vs %s (%s = 0).",
                            label, lab1, col, lab0, col))
}

make_menarche_stratum <- function() {
  e <- check_01(opt$early_menarche_col)
  l <- check_01(opt$late_menarche_col)
  both <- which(!is.na(e) & !is.na(l) & e == 1 & l == 1)
  if (length(both) > 0) {
    stop("Participants flagged as both early and late menarche: ",
         paste(utils::head(ids[both], 10), collapse = ", "))
  }
  g <- rep(NA_character_, length(e))
  g[!is.na(e) & e == 1] <- menarche_labels[1]
  g[!is.na(l) & l == 1] <- menarche_labels[3]
  g[!is.na(e) & !is.na(l) & e == 0 & l == 0] <- menarche_labels[2]
  list(key = "menarche", label = "Menarche timing", sheet = "Menarche",
       group = factor(g, levels = menarche_labels),
       defining = c(opt$early_menarche_col, opt$late_menarche_col,
                    opt$menarche_age_col),
       definition = sprintf(
         "Menarche timing: %s = 1 → %s; %s = 1 → %s; both 0 → %s",
         opt$early_menarche_col, menarche_labels[1],
         opt$late_menarche_col, menarche_labels[3], menarche_labels[2]))
}

strata <- lapply(strata_keys, function(k) {
  switch(k,
    race     = make_binary_stratum("race", "Race", "Race", opt$race_col,
                                   "Non-Black", "Black"),
    obesity  = make_binary_stratum("obesity", "Obesity", "Obesity", opt$obesity_col,
                                   "Not obese", "Obese"),
    menarche = make_menarche_stratum())
})

for (s in strata) {
  tab <- table(s$group, useNA = "always")
  names(tab)[is.na(names(tab))] <- "Missing"
  message("Stratum ", s$label, ": ",
          paste0(names(tab), " = ", as.integer(tab), collapse = "; "))
}

# ── Verification report ────────────────────────────────────────────────────
if (isTRUE(opt$verify)) {
  message("\n", strrep("=", 72), "\nVERIFICATION REPORT\n", strrep("=", 72))
  message(sprintf("Samples: %d   Traits: %d   Strata: %s", n_total, length(traits),
                  paste(vapply(strata, `[[`, "", "label"), collapse = ", ")))
  message(sprintf("\n  %-32s %-11s %6s %6s %10s %10s",
                  "Trait", "Type", "N", "N_NA", "Min", "Max"))
  for (t in traits) {
    x <- dat[[t]]
    message(sprintf("  %-32s %-11s %6d %6d %10.4g %10.4g",
                    substr(t, 1, 32), trait_type[[t]], sum(!is.na(x)),
                    sum(is.na(x)), min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  }
  if (length(strata) > 1) {
    message("\nCross-tabulation of strata (first two):")
    print(table(strata[[1]]$group, strata[[2]]$group, useNA = "ifany",
                dnn = c(strata[[1]]$label, strata[[2]]$label)))
  }
  message(strrep("=", 72), "\n")
}

# ── Formatting and tests ───────────────────────────────────────────────────
DASH <- "—"
NDASH <- "–"

# Continuous summaries are formatted with ONE number of decimals per row, so
# every cell of a trait lines up (52.0 (1.4) / 52.1 (1.5), not 52 / 52.1).
# Decimals come from the trait's overall centre and --digits significant
# digits: age ~52 -> 1 decimal; DDT ~12 -> 1; PCB ~0.07 -> 4 (0.0712).
# Integer-valued traits (counts, scores) get 0 decimals for medians and
# quartiles when those are whole numbers, else 1; 1 decimal for mean (SD).
cont_stats <- function(x, use_mean_sd) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NULL)
  if (use_mean_sd) {
    c(mean(x), if (length(x) > 1) stats::sd(x) else NA_real_)
  } else {
    stats::quantile(x, c(0.5, 0.25, 0.75), names = FALSE, type = 7)
  }
}

choose_decimals <- function(x, stat_list, use_mean_sd) {
  x <- x[!is.na(x)]
  vals <- unlist(stat_list)
  vals <- vals[is.finite(vals)]
  if (all(x == round(x))) {
    if (use_mean_sd) return(1L)
    return(if (all(vals == round(vals))) 0L else 1L)
  }
  ref <- abs(if (use_mean_sd) mean(x) else stats::median(x))
  if (ref == 0) ref <- max(abs(x))
  if (ref == 0) return(0L)
  as.integer(min(6, max(0, opt$digits - 1 - floor(log10(ref)))))
}

fmt_cont_cell <- function(s, d, use_mean_sd) {
  if (is.null(s)) return(DASH)
  f <- function(v) formatC(v, format = "f", digits = d)
  if (use_mean_sd) {
    if (is.na(s[2])) return(f(s[1]))
    sprintf("%s (%s)", f(s[1]), f(s[2]))
  } else {
    sprintf("%s (%s%s%s)", f(s[1]), f(s[2]), NDASH, f(s[3]))
  }
}

summ_bin <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(DASH)
  sprintf("%d (%.1f%%)", sum(x == 1), 100 * mean(x == 1))
}

test_notes <- character(0)

p_cont <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  if (nlevels(g) < 2 || any(table(g) == 0) || length(unique(x)) < 2) return(NA_real_)
  if (nlevels(g) == 2) {
    stats::wilcox.test(x ~ g, exact = FALSE)$p.value
  } else {
    stats::kruskal.test(x ~ g)$p.value
  }
}

p_bin <- function(x, g, trait, stratum) {
  ok <- !is.na(x) & !is.na(g)
  x <- factor(x[ok], levels = c(0, 1)); g <- droplevels(g[ok])
  if (nlevels(g) < 2 || length(unique(x)) < 2) return(NA_real_)
  tab <- table(x, g)
  if (opt$binary_test == "chisq") {
    return(suppressWarnings(stats::chisq.test(tab, correct = FALSE)$p.value))
  }
  res <- tryCatch(stats::fisher.test(tab)$p.value, error = function(e) e)
  if (inherits(res, "error")) {
    # Exact network algorithm can exhaust its workspace on larger r x c
    # tables; fall back to a seeded Monte Carlo exact p and record it.
    set.seed(20261021)
    res <- stats::fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
    test_notes <<- c(test_notes, sprintf(
      "%s x %s: Fisher exact p estimated by Monte Carlo (B = 1e5, seed 20261021).",
      trait, stratum))
  }
  res
}

fmt_p <- function(p) {
  if (is.na(p)) return(DASH)
  if (p < 0.001) return("<0.001")
  if (p < 0.01)  return(sprintf("%.3f", p))
  sprintf("%.2f", p)
}

# ── Build the table body ───────────────────────────────────────────────────
# Each row: list(kind = "section"|"trait", label, cells (character), sig (logical))
# Column layout: Variable | N | Overall | for each stratum: levels..., p
level_counts <- lapply(strata, function(s) as.integer(table(s$group)))

col_spec <- data.frame(key = c("variable", "n", "overall"),
                       stratum = c(NA, NA, NA),
                       header = c("Characteristic", "N",
                                  sprintf("Overall\n(N = %d)", n_total)),
                       is_p = FALSE, stringsAsFactors = FALSE)
for (i in seq_along(strata)) {
  s <- strata[[i]]
  col_spec <- rbind(col_spec, data.frame(
    key = paste0(s$key, "_", seq_along(levels(s$group))),
    stratum = s$key,
    header = sprintf("%s\n(N = %d)", levels(s$group), level_counts[[i]]),
    is_p = FALSE, stringsAsFactors = FALSE))
  col_spec <- rbind(col_spec, data.frame(key = paste0(s$key, "_p"), stratum = s$key,
                                         header = "p", is_p = TRUE,
                                         stringsAsFactors = FALSE))
}

rows <- list()
pvals_long <- list()
prev_group <- NULL
for (i in seq_along(traits)) {
  t <- traits[i]
  g_i <- groups[i]
  if (!is.na(g_i) && !identical(g_i, prev_group)) {
    rows[[length(rows) + 1]] <- list(kind = "section", label = g_i)
    prev_group <- g_i
  }
  x    <- dat[[t]]
  type <- trait_type[[t]]
  use_msd <- t %in% mean_sd_traits || (type == "continuous" && opt$cont_stat == "mean_sd")
  if (type == "binary") {
    summ <- summ_bin
  } else {
    # Decide decimals once per row from every summary that will be shown
    subsets <- c(list(x), unlist(lapply(strata, function(s) {
      lapply(levels(s$group), function(l) x[!is.na(s$group) & s$group == l])
    }), recursive = FALSE))
    d_row <- choose_decimals(x, lapply(subsets, cont_stats, use_mean_sd = use_msd),
                             use_msd)
    summ <- function(v) fmt_cont_cell(cont_stats(v, use_msd), d_row, use_msd)
  }
  stat_lbl <- if (type == "binary") "n (%)" else if (use_msd) "mean (SD)" else
    paste0("median (Q1", NDASH, "Q3)")

  cells <- c(t, as.character(sum(!is.na(x))), summ(x))
  sig   <- c(FALSE, FALSE, FALSE)
  for (s in strata) {
    lv <- levels(s$group)
    if (t %in% s$defining) {
      cells <- c(cells, rep(DASH, length(lv) + 1))
      sig   <- c(sig, rep(FALSE, length(lv) + 1))
      pvals_long[[length(pvals_long) + 1]] <- data.frame(
        Trait = t, Stratum = s$label, Test = "not tested (defining variable)",
        P = NA_real_, stringsAsFactors = FALSE)
      next
    }
    cells <- c(cells, vapply(lv, function(l) summ(x[!is.na(s$group) & s$group == l]),
                             character(1)))
    p <- if (type == "binary") p_bin(x, s$group, t, s$label) else p_cont(x, s$group)
    cells <- c(cells, fmt_p(p))
    sig   <- c(sig, rep(FALSE, length(lv)), isTRUE(opt$highlight_sig) && !is.na(p) && p < opt$alpha)
    test_name <- if (type == "binary") {
      if (opt$binary_test == "fisher") "Fisher exact" else "Pearson chi-squared"
    } else if (length(lv) == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis"
    pvals_long[[length(pvals_long) + 1]] <- data.frame(
      Trait = t, Stratum = s$label, Test = test_name, P = p,
      stringsAsFactors = FALSE)
  }
  rows[[length(rows) + 1]] <- list(kind = "trait", label = t, cells = cells,
                                   sig = sig, type = type, stat = stat_lbl)
}
pvals_long <- do.call(rbind, pvals_long)

# Row label carries its summary statistic, e.g. "DDT, median (Q1–Q3)"
for (i in seq_along(rows)) {
  if (rows[[i]]$kind == "trait") {
    rows[[i]]$cells[1] <- paste0(rows[[i]]$label, ", ", rows[[i]]$stat)
  }
}

# ── Footnotes ──────────────────────────────────────────────────────────────
cont_note <- if (opt$cont_stat == "median_iqr") {
  paste0("Continuous variables: median (Q1", NDASH, "Q3)",
         if (length(mean_sd_traits)) " unless marked mean (SD)" else "", ".")
} else "Continuous variables: mean (SD)."
footnotes <- c(
  paste(cont_note, "Binary variables: n (%) of participants with non-missing data;",
        "N = participants with non-missing data."),
  paste0("p-values: Wilcoxon rank-sum (2 groups) or Kruskal", NDASH,
         "Wallis (3 groups) for continuous variables; ",
         if (opt$binary_test == "fisher") "Fisher's exact test" else
           "Pearson's chi-squared test", " for binary variables.",
         if (isTRUE(opt$highlight_sig)) sprintf(" Bold: p < %g.", opt$alpha) else ""),
  paste0(DASH, " in a stratum's columns: variable defines that stratum and is not tested."),
  vapply(strata, `[[`, "", "definition"),
  "ICE, Index of Concentration at the Extremes; Q1 and Q4, lowest and highest quartiles.",
  if (length(recoded)) paste0("Values coded ", paste(shQuote(na_codes, "sh"), collapse = "/"),
                              " (failed assay) treated as missing."),
  test_notes
)

# ── Assemble tables ────────────────────────────────────────────────────────
build_table <- function(stratum_keys) {
  keep <- col_spec$key %in% c("variable", "n", "overall") |
    col_spec$stratum %in% stratum_keys
  spec <- col_spec[keep, ]
  body <- do.call(rbind, lapply(rows, function(r) {
    if (r$kind == "section") c(r$label, rep("", sum(keep) - 1)) else r$cells[keep]
  }))
  sig <- do.call(rbind, lapply(rows, function(r) {
    if (r$kind == "section") rep(FALSE, sum(keep)) else r$sig[keep]
  }))
  kind <- vapply(rows, `[[`, "", "kind")
  list(spec = spec, body = body, sig = sig, kind = kind)
}

tables <- c(list(build_table(strata_keys)),
            lapply(strata_keys, build_table))
table_names  <- c("Table1", vapply(strata, `[[`, "", "sheet"))
table_titles <- c(opt$title,
                  paste0(opt$title, ", by ", tolower(vapply(strata, `[[`, "", "label"))))

stratum_label <- stats::setNames(vapply(strata, `[[`, "", "label"), strata_keys)

# ── Excel workbook ─────────────────────────────────────────────────────────
xlsx_path <- file.path(opt$out_dir, paste0(opt$prefix, ".xlsx"))
message("Writing workbook: ", xlsx_path)

wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 10, fontName = "Arial")
st_title   <- createStyle(textDecoration = "bold", fontSize = 11)
st_span    <- createStyle(textDecoration = "bold", halign = "center",
                          border = "bottom", borderStyle = "thin")
st_head    <- createStyle(textDecoration = "bold", halign = "center",
                          valign = "center", wrapText = TRUE,
                          border = "bottom", borderStyle = "thin")
st_head_l  <- createStyle(textDecoration = "bold", halign = "left",
                          valign = "center", wrapText = TRUE,
                          border = "bottom", borderStyle = "thin")
st_top     <- createStyle(border = "top", borderStyle = "thin")
st_bottom  <- createStyle(border = "bottom", borderStyle = "thin")
st_section <- createStyle(textDecoration = c("bold", "italic"), halign = "left")
st_label   <- createStyle(halign = "left", indent = 1)
st_cell    <- createStyle(halign = "center")
st_sig     <- createStyle(textDecoration = "bold")
st_note    <- createStyle(fontSize = 8, textDecoration = "italic")

write_table_sheet <- function(sheet, tb, title) {
  addWorksheet(wb, sheet, gridLines = FALSE)
  spec <- tb$spec; body <- tb$body; nc <- ncol(body); nr <- nrow(body)

  writeData(wb, sheet, title, startRow = 1, startCol = 1)
  addStyle(wb, sheet, st_title, rows = 1, cols = 1)

  # Row 2: stratum spanners; row 3: column headers
  addStyle(wb, sheet, st_top, rows = 2, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
  for (k in unique(stats::na.omit(spec$stratum))) {
    cols <- which(spec$stratum == k)
    writeData(wb, sheet, stratum_label[[k]], startRow = 2, startCol = min(cols))
    mergeCells(wb, sheet, rows = 2, cols = cols)
    addStyle(wb, sheet, st_span, rows = 2, cols = cols, gridExpand = TRUE, stack = TRUE)
  }
  writeData(wb, sheet, as.data.frame(t(spec$header)), startRow = 3,
            colNames = FALSE)
  addStyle(wb, sheet, st_head, rows = 3, cols = 2:nc, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet, st_head_l, rows = 3, cols = 1, stack = TRUE)
  setRowHeights(wb, sheet, rows = 3, heights = 30)

  # Body
  first <- 4
  writeData(wb, sheet, as.data.frame(body, stringsAsFactors = FALSE),
            startRow = first, colNames = FALSE)
  body_rows <- first:(first + nr - 1)
  addStyle(wb, sheet, st_cell, rows = body_rows, cols = 2:nc, gridExpand = TRUE, stack = TRUE)
  for (i in seq_len(nr)) {
    r <- first + i - 1
    if (tb$kind[i] == "section") {
      mergeCells(wb, sheet, rows = r, cols = 1:nc)
      addStyle(wb, sheet, st_section, rows = r, cols = 1, stack = TRUE)
    } else {
      addStyle(wb, sheet, st_label, rows = r, cols = 1, stack = TRUE)
      sc <- which(tb$sig[i, ])
      if (length(sc)) addStyle(wb, sheet, st_sig, rows = r, cols = sc, stack = TRUE)
    }
  }
  addStyle(wb, sheet, st_bottom, rows = max(body_rows), cols = 1:nc,
           gridExpand = TRUE, stack = TRUE)

  # Footnotes
  fn_rows <- max(body_rows) + 1 + seq_along(footnotes)
  writeData(wb, sheet, footnotes, startRow = min(fn_rows), startCol = 1)
  addStyle(wb, sheet, st_note, rows = fn_rows, cols = 1, stack = TRUE)

  setColWidths(wb, sheet, cols = 1, widths = 42)
  setColWidths(wb, sheet, cols = 2, widths = 6)
  setColWidths(wb, sheet, cols = which(!spec$is_p)[-(1:2)], widths = 17)
  setColWidths(wb, sheet, cols = which(spec$is_p), widths = 8)
  freezePane(wb, sheet, firstActiveRow = first, firstActiveCol = 2)
}

for (i in seq_along(tables)) write_table_sheet(table_names[i], tables[[i]], table_titles[i])

# Variable_info
var_info <- data.frame(
  Trait   = traits,
  Section = groups,
  Type    = unname(trait_type[traits]),
  Summary = vapply(rows[vapply(rows, `[[`, "", "kind") == "trait"], `[[`, "", "stat"),
  N_nonmissing = vapply(traits, function(t) sum(!is.na(dat[[t]])), integer(1)),
  N_missing    = vapply(traits, function(t) sum(is.na(dat[[t]])), integer(1)),
  N_recoded_na_codes = unname(co$recoded[traits]),
  N_distinct   = vapply(traits, function(t) length(unique(stats::na.omit(dat[[t]]))), integer(1)),
  check.names = FALSE, stringsAsFactors = FALSE
)
hdr <- createStyle(textDecoration = "bold")
addWorksheet(wb, "Variable_info")
writeData(wb, "Variable_info", var_info, headerStyle = hdr)
setColWidths(wb, "Variable_info", cols = seq_len(ncol(var_info)), widths = "auto")
freezePane(wb, "Variable_info", firstRow = TRUE)

addWorksheet(wb, "P_values")
pv_out <- pvals_long
pv_out$P_formatted <- vapply(pv_out$P, fmt_p, character(1))
writeData(wb, "P_values", pv_out, headerStyle = hdr)
addStyle(wb, "P_values", createStyle(numFmt = "0.00E+00"),
         rows = seq_len(nrow(pv_out)) + 1, cols = 4)
setColWidths(wb, "P_values", cols = seq_len(ncol(pv_out)), widths = "auto")
freezePane(wb, "P_values", firstRow = TRUE)

strata_counts <- vapply(strata, function(s) {
  tab <- table(s$group, useNA = "always")
  names(tab)[is.na(names(tab))] <- "missing"
  paste0(s$label, ": ", paste0(names(tab), " ", as.integer(tab), collapse = ", "))
}, character(1))

run_info <- data.frame(
  Parameter = c("Input", "Sheet", "ID column", "Trait file", "Samples",
                "Traits", "Binary traits", "Continuous traits",
                "Strata", "Strata counts", "Continuous summary",
                "Mean (SD) traits", "Binary test", "Alpha", "NA codes",
                "Cells recoded to NA", "Run at"),
  Value = c(normalizePath(opt$input), opt$sheet, id_col, normalizePath(opt$traits),
            n_total, length(traits), sum(trait_type == "binary"),
            sum(trait_type == "continuous"),
            paste(vapply(strata, `[[`, "", "label"), collapse = ", "),
            paste(strata_counts, collapse = " | "), opt$cont_stat,
            if (length(mean_sd_traits)) paste(mean_sd_traits, collapse = "; ") else "none",
            opt$binary_test, opt$alpha,
            if (length(na_codes)) paste(na_codes, collapse = ", ") else "none",
            if (length(recoded)) paste0(names(recoded), " (", recoded, ")", collapse = "; ") else "none",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  stringsAsFactors = FALSE)
addWorksheet(wb, "Run_info")
writeData(wb, "Run_info", run_info, headerStyle = hdr)
setColWidths(wb, "Run_info", cols = 1:2, widths = c(24, 100))

saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# ── PDF ────────────────────────────────────────────────────────────────────
font_family <- setup_publication_font(candidates = c("Arial", "Helvetica"),
                                      strict = TRUE)

make_table_grob <- function(tb, title) {
  body <- tb$body; spec <- tb$spec; nr <- nrow(body); nc <- ncol(body)
  df <- as.data.frame(body, stringsAsFactors = FALSE)
  names(df) <- spec$header
  is_sec <- tb$kind == "section"

  # Cell-wise parameters are recycled column-major
  face <- matrix("plain", nr, nc)
  face[is_sec, 1] <- "bold.italic"
  face[tb$sig] <- "bold"
  hj <- matrix(0.5, nr, nc); hj[, 1] <- 0
  xx <- matrix(0.5, nr, nc); xx[, 1] <- 0.02
  xx[!is_sec, 1] <- 0.06   # indent trait rows under section headers

  th <- ttheme_minimal(
    base_family = font_family,
    core = list(fg_params = list(fontsize = 7.5, fontface = as.vector(face),
                                 hjust = as.vector(hj), x = as.vector(xx),
                                 fontfamily = font_family),
                bg_params = list(fill = NA),
                padding = unit(c(3, 1.6), "mm")),
    colhead = list(fg_params = list(fontsize = 7.5, fontface = "bold",
                                    hjust = c(0, rep(0.5, nc - 1)),
                                    x = c(0.02, rep(0.5, nc - 1)),
                                    fontfamily = font_family),
                   bg_params = list(fill = NA),
                   padding = unit(c(3, 2), "mm")))
  g <- tableGrob(df, rows = NULL, theme = th)

  rule <- function(y) segmentsGrob(x0 = 0, x1 = 1, y0 = y, y1 = y,
                                   gp = gpar(lwd = 0.8))
  # Spanner row above the column headers
  g <- gtable_add_rows(g, unit(5, "mm"), pos = 0)
  for (k in unique(stats::na.omit(spec$stratum))) {
    cols <- which(spec$stratum == k)
    g <- gtable_add_grob(g, textGrob(stratum_label[[k]],
                                     gp = gpar(fontsize = 7.5, fontface = "bold",
                                               fontfamily = font_family)),
                         t = 1, l = min(cols), r = max(cols), clip = "off",
                         name = paste0("span_", k))
    g <- gtable_add_grob(g, segmentsGrob(x0 = unit(1.5, "mm"),
                                         x1 = unit(1, "npc") - unit(1.5, "mm"),
                                         y0 = 0, y1 = 0, gp = gpar(lwd = 0.5)),
                         t = 1, l = min(cols), r = max(cols),
                         name = paste0("span_rule_", k))
  }
  g <- gtable_add_grob(g, rule(1), t = 1, l = 1, r = nc, name = "rule_top")
  g <- gtable_add_grob(g, rule(0), t = 2, l = 1, r = nc, name = "rule_head")
  g <- gtable_add_grob(g, rule(0), t = nrow(g), l = 1, r = nc, name = "rule_bottom")

  # Title above, footnotes below
  g <- gtable_add_rows(g, unit(8, "mm"), pos = 0)
  g <- gtable_add_grob(g, textGrob(title, x = 0, hjust = 0,
                                   gp = gpar(fontsize = 10, fontface = "bold",
                                             fontfamily = font_family)),
                       t = 1, l = 1, r = nc, clip = "off", name = "title")

  width_in  <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
  chars_per_line <- max(60, floor(width_in * 72 / (6.5 * 0.5)))
  wrapped <- unlist(lapply(footnotes, strwrap, width = chars_per_line))
  g <- gtable_add_rows(g, unit(2, "mm"))
  for (ln in wrapped) {
    g <- gtable_add_rows(g, unit(3.4, "mm"))
    g <- gtable_add_grob(g, textGrob(ln, x = 0, hjust = 0,
                                     gp = gpar(fontsize = 6.5, fontface = "italic",
                                               fontfamily = font_family)),
                         t = nrow(g), l = 1, r = nc, clip = "off")
  }
  g
}

pdf_path <- file.path(opt$out_dir, paste0(opt$prefix, ".pdf"))
# Build and measure the grobs on a cairo device: the default pdf() device
# cannot encode the em/en dashes and would mis-measure text widths.
measure_file <- tempfile(fileext = ".pdf")
cairo_pdf(measure_file, width = 30, height = 30, family = font_family)
grobs <- Map(make_table_grob, tables, table_titles)
sizes <- lapply(grobs, function(g) c(
  w = convertWidth(sum(g$widths), "in", valueOnly = TRUE),
  h = convertHeight(sum(g$heights), "in", valueOnly = TRUE)))
page_w <- max(vapply(sizes, `[[`, 0, "w")) + 1
page_h <- max(vapply(sizes, `[[`, 0, "h")) + 1
invisible(dev.off())
unlink(measure_file)

cairo_pdf(pdf_path, width = page_w, height = page_h, family = font_family,
          onefile = TRUE)
for (g in grobs) {
  grid.newpage()
  pushViewport(viewport(x = unit(0.5, "in"), y = unit(1, "npc") - unit(0.5, "in"),
                        width = sum(g$widths), height = sum(g$heights),
                        just = c("left", "top")))
  grid.draw(g)
  popViewport()
}
invisible(dev.off())
message("Saved: ", pdf_path, sprintf(" (%.1f x %.1f in)", page_w, page_h))

# ── run_parameters.txt ─────────────────────────────────────────────────────
opt_log <- opt[setdiff(names(opt), "help")]
log_lines <- c(
  "01_table1.R — run parameters",
  paste("Run at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Script dir:", script_dir),
  "",
  "## Arguments",
  sprintf("%s = %s", names(opt_log),
          vapply(opt_log, function(v) paste(v, collapse = ","), character(1))),
  "",
  "## Data",
  paste("Samples:", n_total),
  paste("Traits:", length(traits), sprintf("(%d binary, %d continuous)",
                                           sum(trait_type == "binary"),
                                           sum(trait_type == "continuous"))),
  paste("Cells recoded to NA:",
        if (length(recoded)) paste0(names(recoded), " (", recoded, ")", collapse = "; ") else "none"),
  strata_counts,
  if (length(test_notes)) c("", "## Test notes", test_notes),
  "",
  "## Software",
  R.version.string,
  paste0("openxlsx ", utils::packageVersion("openxlsx"),
         "; gridExtra ", utils::packageVersion("gridExtra"),
         "; gtable ", utils::packageVersion("gtable"))
)
writeLines(log_lines, file.path(opt$out_dir, "run_parameters.txt"))

message("Done. ", length(traits), " traits x ", length(strata), " strata.")