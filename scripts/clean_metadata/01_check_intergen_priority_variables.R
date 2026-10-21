#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Convert selected Intergen metadata variables to numeric
# Keeps full dataset, rownames, and conversion report
# ============================================================

# -----------------------------
# Input / output paths
# -----------------------------
infile <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_070702026_publication_labels.xlsx"

outdir <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/numeric_conversion_check"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile_full_numeric <- file.path(outdir, "Intergen173_070702026_full_dataset_priority_numeric.xlsx")
outfile_report       <- file.path(outdir, "Intergen173_070702026_priority_variable_conversion_report.xlsx")
outfile_missing      <- file.path(outdir, "Intergen173_070702026_priority_variables_missing.txt")


# infile <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_publication_labels.xlsx"

# outdir <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/numeric_conversion_check"
# dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# outfile_full_numeric <- file.path(outdir, "Intergen173_05222026_full_dataset_priority_numeric.xlsx")
# outfile_report       <- file.path(outdir, "Intergen173_05222026_priority_variable_conversion_report.xlsx")
# outfile_missing      <- file.path(outdir, "Intergen173_05222026_priority_variables_missing.txt")

# -----------------------------
# Priority variables to convert/check
# -----------------------------
priority_vars <- unique(c(
  "PCB82", "PCB174", "PCB177", "PCB178", "PCB199",
  "PFHxA", "DDT", "DDE", "op_DDT", "MePFOSAAcOH",
  "PCB201", "PFDeA", "F1 bBHC", "F1 APOE genotype",
  "F1 sex", "PFDoA", "PFNA"
))

# Keep genotype categorical unless manually recoded later
categorical_vars <- c("F1 APOE genotype")

# -----------------------------
# Read Excel and use first column as rownames
# -----------------------------
dat_raw <- readxl::read_excel(infile)
dat_raw <- as.data.frame(dat_raw, check.names = FALSE)

id_col <- names(dat_raw)[1]

if (any(is.na(dat_raw[[id_col]]) | dat_raw[[id_col]] == "")) {
  stop("The first column used for rownames contains missing or blank IDs: ", id_col)
}

rownames(dat_raw) <- make.unique(as.character(dat_raw[[id_col]]))

# Remove ID column from working dataset
dat <- dat_raw[, setdiff(names(dat_raw), id_col), drop = FALSE]

cat("Loaded:", nrow(dat), "rows x", ncol(dat), "columns\n")
cat("Used rownames from first column:", id_col, "\n")

# -----------------------------
# Check priority variables
# -----------------------------
found_vars <- intersect(priority_vars, names(dat))
missing_vars <- setdiff(priority_vars, names(dat))

cat("\nPriority variables requested:", length(priority_vars), "\n")
cat("Found:", length(found_vars), "\n")
cat("Missing:", length(missing_vars), "\n")

if (length(missing_vars) > 0) {
  cat("\nMissing variables:\n")
  print(missing_vars)
}

writeLines(missing_vars, outfile_missing)

# -----------------------------
# Safe numeric conversion function
# -----------------------------
safe_numeric_convert <- function(x) {

  if (is.numeric(x)) {
    return(x)
  }

  x_chr <- as.character(x)
  x_chr <- stringr::str_trim(x_chr)

  # Standardize missing values
  x_chr[x_chr %in% c(
    "", "NA", "N/A", "na", "n/a", ".",
    "Missing", "missing", "Unknown", "unknown",
    "ND", "nd", "Non-detect", "non-detect",
    "Non Detect", "non detect"
  )] <- NA

  x_lower <- stringr::str_to_lower(x_chr)

  # Convert common binary labels
  x_mapped <- dplyr::case_when(
    x_lower %in% c("yes", "y", "true", "t", "present", "positive", "male", "m") ~ "1",
    x_lower %in% c("no", "n", "false", "f", "absent", "negative", "female", "fem") ~ "0",
    TRUE ~ x_chr
  )

  suppressWarnings(as.numeric(x_mapped))
}

# -----------------------------
# Create full output dataset
# -----------------------------
dat_clean <- dat

# Only convert the priority variables.
# All other variables remain unchanged.
for (v in found_vars) {

  if (v %in% categorical_vars) {
    dat_clean[[v]] <- as.character(dat[[v]])
  } else {
    dat_clean[[v]] <- safe_numeric_convert(dat[[v]])
  }
}

# Preserve rownames
rownames(dat_clean) <- rownames(dat)

# -----------------------------
# Build conversion report for priority variables only
# -----------------------------
conversion_report <- lapply(found_vars, function(v) {

  x <- dat[[v]]
  x_clean <- dat_clean[[v]]

  non_missing_original <- sum(!is.na(x))
  non_missing_clean <- sum(!is.na(x_clean))

  values_lost <- non_missing_original - non_missing_clean
  prop_lost <- ifelse(non_missing_original == 0, NA, values_lost / non_missing_original)

  data.frame(
    variable = v,
    original_class = paste(class(x), collapse = ";"),
    cleaned_class = paste(class(x_clean), collapse = ";"),
    already_numeric = is.numeric(x),
    kept_categorical = v %in% categorical_vars,
    non_missing_original = non_missing_original,
    non_missing_after_cleaning = non_missing_clean,
    values_lost_during_conversion = values_lost,
    proportion_lost = prop_lost,
    unique_original_values = paste(head(unique(as.character(x)), 30), collapse = " | "),
    unique_cleaned_values = paste(head(unique(as.character(x_clean)), 30), collapse = " | "),
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows()

# -----------------------------
# Quick checks
# -----------------------------
if ("APOE_Genotype" %in% names(dat)) {
  cat("\nAPOE_Genotype original distribution:\n")
  print(table(dat$APOE_Genotype, useNA = "ifany"))

  cat("\nAPOE_Genotype in cleaned dataset, kept categorical:\n")
  print(table(dat_clean$APOE_Genotype, useNA = "ifany"))
}

if ("F1sex" %in% names(dat)) {
  cat("\nF1sex original distribution:\n")
  print(table(dat$F1sex, useNA = "ifany"))

  cat("\nF1sex after numeric conversion:\n")
  print(table(dat_clean$F1sex, useNA = "ifany"))
}

# -----------------------------
# Save full dataset with rownames as first column
# -----------------------------
dat_clean_out <- data.frame(
  rowname = rownames(dat_clean),
  dat_clean,
  check.names = FALSE
)

report_out <- data.frame(
  rowname = seq_len(nrow(conversion_report)),
  conversion_report,
  check.names = FALSE
)

openxlsx::write.xlsx(
  dat_clean_out,
  outfile_full_numeric,
  rowNames = FALSE,
  overwrite = TRUE
)

openxlsx::write.xlsx(
  report_out,
  outfile_report,
  rowNames = FALSE,
  overwrite = TRUE
)

cat("\nSaved full dataset with priority variables converted to numeric:\n")
cat(outfile_full_numeric, "\n")

cat("\nSaved priority-variable conversion report:\n")
cat(outfile_report, "\n")

cat("\nSaved missing variable list:\n")
cat(outfile_missing, "\n")

cat("\nDone.\n")