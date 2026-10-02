# =============================================================================
# Table 1 — Intergen Sample Characteristics
# Stratified by: Race, Obesity, Breast Cancer (F1), Menarche (young/old)
#
# OUTPUT: table1_intergen.xlsx (presentation-ready)
#
# USAGE:
#   Rscript table1_intergen.R --input /path/to/your/data.xlsx
#   OR source interactively and set input_file below
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(gridExtra)
  library(grid)
  library(ggplot2)
})

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Path to input Excel (.xlsx) file [required]"),
  make_option("--sheet", type = "integer", default = 1,
              help = "Sheet number to read [default = 1]"),
  make_option("--output", type = "character", default = "table1_intergen.xlsx",
              help = "Output Excel file path [default = table1_intergen.xlsx]"),
  make_option("--output_pdf", type = "character", default = "table1_intergen.pdf",
              help = "Output PDF file path [default = table1_intergen.pdf]"),

  make_option("--verify", action = "store_true", default = FALSE,
              help = paste("Print a verification report of the input data before",
                           "building the table: column presence, N per group,",
                           "value distributions, and missing counts. [default = FALSE]"))
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) stop("--input is required. Provide path to your .xlsx file.")

# ── Load data ─────────────────────────────────────────────────────────────────
message("Reading: ", opt$input)
df <- read_excel(opt$input, sheet = opt$sheet)
message("Loaded ", nrow(df), " rows, ", ncol(df), " columns")


# ── Variable definitions ──────────────────────────────────────────────────────

# Grouping variables (each creates one stratified block in the table)

# Derive menarche group: 1=young, 2=old, 0=middle (neither flagged)
df <- df %>%
  mutate(
    menarche_group = case_when(
      youngmen_adol == 1 ~ "Young (<12 yrs)",
      oldmen_adol   == 1 ~ "Old (>14 yrs)",
      TRUE               ~ "Middle"
    )
  )

# Row variables: exposures + sample characteristics
# Format: list(display_label = column_name)
continuous_vars <- list(
  "Age at WGBS (years)"              = "Cati_Age",
  "DDE (ng/mL)"                      = "dde",
  "DDT (ng/mL)"                      = "ddt",
  "op-DDT (ng/mL)"                   = "op_ddt",
  "ICE Education Q1 (most disadvantaged)"   = "ICE60edu1_q1",
  "ICE Race Q1"                      = "ICE60rac1_q1",
  "ICE Income Q1"                    = "ICE60inc1_q1",
  "ICE Inc×Race Q1"                  = "ICE60incxrac1_q1",
  "Age at menarche (adolescent)"     = "agemenarche_adol"
)

binary_vars <- list(
  "DDE tertile 1 (lowest)"   = "ddeT1",
  "DDE tertile 2"             = "ddeT2",
  "DDE tertile 3 (highest)"  = "ddeT3",
  "DDT tertile 1 (lowest)"   = "ddtT1",
  "DDT tertile 2"             = "ddtT2",
  "DDT tertile 3 (highest)"  = "ddtT3",
  "Obese F1"                  = "HVobeseF1",
  "Overweight F1"             = "HVoverwtF1",
  "Breast cancer F0"          = "breastcaF0",
  "Breast cancer F1"          = "breastcaF1",
  "Young menarche (<12)"      = "youngmen_adol",
  "Old menarche (>14)"        = "oldmen_adol"
)


# ── Verification report (runs after variable definitions) ─────────────────────
if (isTRUE(opt$verify)) {
  message("\n", strrep("=", 60))
  message("VERIFICATION REPORT")
  message(strrep("=", 60))

  # 1. Column presence check
  expected_cols <- c(
    unlist(continuous_vars), unlist(binary_vars),
    "aakid", "HVobeseF1", "breastcaF1", "youngmen_adol", "oldmen_adol"
  )
  message("\n── Column presence check ──────────────────────────────")
  for (col in unique(expected_cols)) {
    status <- if (col %in% names(df)) "OK  " else "MISSING"
    message(sprintf("  %-35s %s", col, status))
  }

  # 2. Grouping variable distributions
  message("\n── Grouping variable distributions ────────────────────")
  grp_checks <- list(
    "Race      (aakid)"         = "aakid",
    "Obesity   (HVobeseF1)"     = "HVobeseF1",
    "BreastCa  (breastcaF1)"    = "breastcaF1",
    "Menarche  (youngmen_adol)" = "youngmen_adol",
    "Menarche  (oldmen_adol)"   = "oldmen_adol"
  )
  for (lbl in names(grp_checks)) {
    col <- grp_checks[[lbl]]
    if (!col %in% names(df)) next
    tbl <- table(df[[col]], useNA = "always")
    message(sprintf("\n  %s", lbl))
    for (i in seq_along(tbl)) {
      nm      <- names(tbl)[i]
      val_lbl <- if (is.na(nm) || nm == "NA") "NA" else nm
      message(sprintf("    %-10s : %d (%.1f%%)",
                      val_lbl, tbl[i], 100 * tbl[i] / nrow(df)))
    }
  }

  # 3. Continuous variable summaries
  message("\n── Continuous variable summaries ───────────────────────")
  message(sprintf("  %-35s %8s %8s %8s %8s %6s",
                  "Variable", "Min", "Mean", "Max", "SD", "N_NA"))
  for (lbl in names(continuous_vars)) {
    col <- continuous_vars[[lbl]]
    if (!col %in% names(df)) next
    x    <- as.numeric(df[[col]])
    n_na <- sum(is.na(x))
    x_ok <- x[!is.na(x)]
    message(sprintf("  %-35s %8.2f %8.2f %8.2f %8.2f %6d",
                    lbl, min(x_ok), mean(x_ok), max(x_ok), sd(x_ok), n_na))
  }

  # 4. Binary variable summaries
  message("\n── Binary variable summaries (N=1 / total non-NA) ─────")
  message(sprintf("  %-35s %8s %8s %6s", "Variable", "N=1", "Pct=1", "N_NA"))
  for (lbl in names(binary_vars)) {
    col <- binary_vars[[lbl]]
    if (!col %in% names(df)) next
    x    <- as.numeric(df[[col]])
    n_na <- sum(is.na(x))
    n1   <- sum(x == 1, na.rm = TRUE)
    tot  <- sum(!is.na(x))
    message(sprintf("  %-35s %8d %7.1f%% %6d",
                    lbl, n1, 100 * n1 / max(tot, 1), n_na))
  }

  # 5. Derived menarche group
  df_tmp <- df %>%
    mutate(menarche_group = case_when(
      youngmen_adol == 1 ~ "Young (<12 yrs)",
      oldmen_adol   == 1 ~ "Old (>14 yrs)",
      TRUE               ~ "Middle"
    ))
  message("\n── Derived menarche_group ──────────────────────────────")
  tbl_men <- table(df_tmp$menarche_group, useNA = "always")
  for (i in seq_along(tbl_men)) {
    nm      <- names(tbl_men)[i]
    val_lbl <- if (is.na(nm) || nm == "NA") "NA" else nm
    message(sprintf("  %-20s : %d (%.1f%%)",
                    val_lbl, tbl_men[i], 100 * tbl_men[i] / nrow(df)))
  }

  # 6. Cross-check key values
  message("\n── Cross-check: raw data values ────────────────────────")
  checks <- list(
    list(label = "N total",                       val = nrow(df)),
    list(label = "N Black (aakid=1)",              val = sum(df$aakid == 1, na.rm = TRUE)),
    list(label = "N Non-Black (aakid=0)",          val = sum(df$aakid == 0, na.rm = TRUE)),
    list(label = "N Obese (HVobeseF1=1)",          val = sum(df$HVobeseF1 == 1, na.rm = TRUE)),
    list(label = "N Not Obese (HVobeseF1=0)",      val = sum(df$HVobeseF1 == 0, na.rm = TRUE)),
    list(label = "N Breast Cancer F1 (breastcaF1=1)",       val = sum(df$breastcaF1 == 1, na.rm = TRUE)),
    list(label = "N Young menarche",               val = sum(df$youngmen_adol == 1, na.rm = TRUE)),
    list(label = "N Old menarche",                 val = sum(df$oldmen_adol == 1, na.rm = TRUE)),
    list(label = "Mean DDE overall",               val = round(mean(df$dde, na.rm = TRUE), 1)),
    list(label = "Mean DDE Black women",           val = round(mean(df$dde[df$aakid == 1], na.rm = TRUE), 1)),
    list(label = "Mean DDE Non-Black women",       val = round(mean(df$dde[df$aakid == 0], na.rm = TRUE), 1)),
    list(label = "N Obese among Black women",      val = sum(df$HVobeseF1[df$aakid == 1] == 1, na.rm = TRUE)),
    list(label = "N Obese among Non-Black women",  val = sum(df$HVobeseF1[df$aakid == 0] == 1, na.rm = TRUE))
  )
  message(sprintf("  %-40s %12s", "Check", "Value"))
  message(sprintf("  %s", strrep("-", 54)))
  for (chk in checks) {
    message(sprintf("  %-40s %12s", chk$label, as.character(chk$val)))
  }
  message("\n  Verify these match what you see in the output table.")
  message(strrep("=", 60))
  message("END VERIFICATION")
  message(strrep("=", 60), "\n")
}

# ── Helper functions ──────────────────────────────────────────────────────────

fmt_mean_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

fmt_n_pct <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  n   <- sum(x == 1, na.rm = TRUE)
  tot <- length(x)
  if (tot == 0) return("—")
  sprintf("%d (%.1f%%)", n, 100 * n / tot)
}

p_value_cont <- function(x, grp) {
  x   <- as.numeric(x)
  grp <- as.factor(grp)
  valid <- !is.na(x) & !is.na(grp)
  if (sum(valid) < 3 || length(unique(grp[valid])) < 2) return(NA_real_)
  tryCatch({
    if (length(levels(grp[valid, drop = TRUE])) == 2) {
      wilcox.test(x[valid] ~ grp[valid])$p.value
    } else {
      kruskal.test(x[valid] ~ grp[valid])$p.value
    }
  }, error = function(e) NA_real_)
}

p_value_binary <- function(x, grp) {
  x   <- as.numeric(x)
  grp <- as.factor(grp)
  valid <- !is.na(x) & !is.na(grp)
  if (sum(valid) < 3 || length(unique(grp[valid])) < 2) return(NA_real_)
  tryCatch({
    tbl <- table(x[valid], grp[valid])
    chisq.test(tbl, simulate.p.value = TRUE)$p.value
  }, error = function(e) NA_real_)
}

fmt_p <- function(p) {
  if (is.na(p)) return("—")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# ── Build one stratified block per grouping variable ──────────────────────────
# Every table always includes fixed extra columns for:
#   - Black vs Non-Black (race)
#   - Obese vs Not Obese
# so these are visible in every stratification regardless of primary group.

build_block <- function(data, group_col, group_label) {
  grp_vec   <- data[[group_col]]
  n_overall <- nrow(data)

  # Ordered group levels and labels for publication-ready output
  grp_meta <- switch(
    group_col,
    "aakid" = list(levels = c(0, 1), labels = c("Non-Black", "Black")),
    "HVobeseF1" = list(levels = c(0, 1), labels = c("Not Obese", "Obese")),
    "breastcaF1" = list(levels = c(0, 1), labels = c("No BrCa", "BrCa F1")),
    "menarche_group" = list(levels = c("Middle", "Young (<12 yrs)", "Old (>14 yrs)"),
                            labels = c("Middle", "Young (<12 yrs)", "Old (>14 yrs)")),
    list(levels = sort(unique(grp_vec[!is.na(grp_vec)])),
         labels = as.character(sort(unique(grp_vec[!is.na(grp_vec)]))))
  )

  # Keep only levels that occur in the data, but preserve meaningful order
  grp_lvls <- grp_meta$levels[grp_meta$levels %in% unique(grp_vec[!is.na(grp_vec)])]
  grp_display <- grp_meta$labels[grp_meta$levels %in% grp_lvls]
  n_per_grp <- sapply(grp_lvls, function(l) sum(grp_vec == l, na.rm = TRUE))

  col_headers <- c(
    "Variable",
    sprintf("Overall\n(N=%d)", n_overall),
    sprintf("%s\n(N=%d)", grp_display, n_per_grp),
    "p-value"
  )

  make_row <- function(lbl, col, is_binary = FALSE) {
    fmt_fn  <- if (is_binary) fmt_n_pct      else fmt_mean_sd
    pval_fn <- if (is_binary) p_value_binary else p_value_cont
    if (!col %in% names(data)) return(NULL)

    c(paste0("    ", lbl),
      fmt_fn(data[[col]]),
      sapply(grp_lvls, function(l)
        fmt_fn(data[[col]][grp_vec == l & !is.na(grp_vec)])),
      fmt_p(pval_fn(data[[col]], grp_vec)))
  }

  rows <- list()
  rows[[length(rows)+1]] <- c(paste0("── Stratified by: ", group_label),
                               rep("", length(col_headers) - 1))
  rows[[length(rows)+1]] <- c("  N", as.character(n_overall),
                               as.character(n_per_grp), "")
  rows[[length(rows)+1]] <- c("  Continuous — Mean (SD)",
                               rep("", length(col_headers) - 1))
  for (lbl in names(continuous_vars)) {
    r <- make_row(lbl, continuous_vars[[lbl]], FALSE)
    if (!is.null(r)) rows[[length(rows)+1]] <- r
  }
  rows[[length(rows)+1]] <- c("  Binary — N (%)",
                               rep("", length(col_headers) - 1))
  for (lbl in names(binary_vars)) {
    r <- make_row(lbl, binary_vars[[lbl]], TRUE)
    if (!is.null(r)) rows[[length(rows)+1]] <- r
  }

  tbl <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(tbl) <- col_headers
  tbl
}

# ── Build all blocks ──────────────────────────────────────────────────────────

blocks <- list(
  build_block(df, "aakid",        "Race (Black vs Non-Black)"),
  build_block(df, "HVobeseF1",    "Obesity"),
  build_block(df, "breastcaF1",   "Breast Cancer F1"),
  build_block(df, "menarche_group", "Menarche Timing")
)

# ── Write to Excel ────────────────────────────────────────────────────────────

wb <- createWorkbook()

# Styles
style_title   <- createStyle(fontSize = 13, fontColour = "#FFFFFF",
                              fgFill = "#1F4E79", fontName = "Arial",
                              textDecoration = "bold", halign = "center",
                              wrapText = TRUE, valign = "center")
style_section <- createStyle(fontSize = 10, fontColour = "#FFFFFF",
                              fgFill = "#2E75B6", fontName = "Arial",
                              textDecoration = "bold")
style_subhdr  <- createStyle(fontSize = 9,  fontColour = "#1F4E79",
                              fgFill = "#D6E4F0", fontName = "Arial",
                              textDecoration = "bold", halign = "left")
style_colhdr  <- createStyle(fontSize = 9,  fontColour = "#FFFFFF",
                              fgFill = "#2E75B6", fontName = "Arial",
                              textDecoration = "bold", halign = "center",
                              wrapText = TRUE, valign = "center")
style_data    <- createStyle(fontSize = 9,  fontName = "Arial",
                              halign = "center", border = "TopBottomLeftRight",
                              borderColour = "#D9D9D9")
style_rowlbl  <- createStyle(fontSize = 9,  fontName = "Arial",
                              halign = "left",   border = "TopBottomLeftRight",
                              borderColour = "#D9D9D9")
style_alt     <- createStyle(fontSize = 9,  fontName = "Arial",
                              fgFill = "#F2F7FC", halign = "center",
                              border = "TopBottomLeftRight",
                              borderColour = "#D9D9D9")
style_altlbl  <- createStyle(fontSize = 9,  fontName = "Arial",
                              fgFill = "#F2F7FC", halign = "left",
                              border = "TopBottomLeftRight",
                              borderColour = "#D9D9D9")
style_sig     <- createStyle(fontSize = 9,  fontName = "Arial",
                              fontColour = "#8B0000", textDecoration = "bold",
                              halign = "center", border = "TopBottomLeftRight",
                              borderColour = "#D9D9D9")

# ── Sheet 0: Manuscript-ready combined summary ────────────────────────────────
# One compact sheet showing Overall + each primary comparison once.
# This avoids repeating Race/Obesity columns inside every separate stratified sheet.
addWorksheet(wb, "Manuscript_Table1")
ws_all <- "Manuscript_Table1"

summary_specs <- list(
  list(name = "race", label = "race", col = "aakid",
       levels = c(0, 1), labels = c("Non-Black", "Black")),
  list(name = "obesity", label = "obesity", col = "HVobeseF1",
       levels = c(0, 1), labels = c("Not Obese", "Obese")),
  list(name = "Breast Cancer", label = "Breast Cancer F1", col = "breastcaF1",
       levels = c(0, 1), labels = c("No BrCa F1", "BrCa F1")),
  list(name = "menarche", label = "menarche", col = "menarche_group",
       levels = c("Middle", "Young (<12 yrs)", "Old (>14 yrs)"),
       labels = c("Middle", "Young (<12 yrs)", "Old (>14 yrs)"))
)

make_summary_headers <- function(data) {
  headers <- c("Variable", sprintf("Overall\n(N=%d)", nrow(data)))
  for (sp in summary_specs) {
    grp <- data[[sp$col]]
    keep <- sp$levels %in% unique(grp[!is.na(grp)])
    lvls <- sp$levels[keep]
    labs <- sp$labels[keep]
    ns <- sapply(lvls, function(l) sum(grp == l, na.rm = TRUE))
    headers <- c(headers, sprintf("%s\n(N=%d)", labs, ns), sprintf("p (%s)", sp$label))
  }
  headers
}

make_summary_row <- function(lbl, col, is_binary = FALSE) {
  fmt_fn  <- if (is_binary) fmt_n_pct      else fmt_mean_sd
  pval_fn <- if (is_binary) p_value_binary else p_value_cont
  if (!col %in% names(df)) return(NULL)

  vals <- c(paste0("    ", lbl), fmt_fn(df[[col]]))
  for (sp in summary_specs) {
    grp <- df[[sp$col]]
    keep <- sp$levels %in% unique(grp[!is.na(grp)])
    lvls <- sp$levels[keep]
    vals <- c(
      vals,
      sapply(lvls, function(l) fmt_fn(df[[col]][grp == l & !is.na(grp)])),
      fmt_p(pval_fn(df[[col]], grp))
    )
  }
  vals
}

sg_headers <- make_summary_headers(df)
n_sg_cols <- length(sg_headers)

sum_rows <- list()
sum_rows[[1]] <- c("  Continuous — Mean (SD)", rep("", n_sg_cols - 1))
for (lbl in names(continuous_vars)) {
  r <- make_summary_row(lbl, continuous_vars[[lbl]], FALSE)
  if (!is.null(r)) sum_rows[[length(sum_rows)+1]] <- r
}
sum_rows[[length(sum_rows)+1]] <- c("  Binary — N (%)", rep("", n_sg_cols - 1))
for (lbl in names(binary_vars)) {
  r <- make_summary_row(lbl, binary_vars[[lbl]], TRUE)
  if (!is.null(r)) sum_rows[[length(sum_rows)+1]] <- r
}

sum_tbl <- as.data.frame(do.call(rbind, sum_rows), stringsAsFactors = FALSE)
colnames(sum_tbl) <- sg_headers
ncols_all <- ncol(sum_tbl)

# Write title
writeData(wb, ws_all,
          data.frame(x = "Table 1 — Intergen Sample Characteristics"),
          startRow = 1, startCol = 1, colNames = FALSE)
mergeCells(wb, ws_all, cols = 1:ncols_all, rows = 1)
addStyle(wb, ws_all, style_title, rows = 1, cols = 1:ncols_all, gridExpand = TRUE)
setRowHeights(wb, ws_all, rows = 1, heights = 30)

# Column headers
writeData(wb, ws_all, as.data.frame(t(colnames(sum_tbl))),
          startRow = 2, startCol = 1, colNames = FALSE)
addStyle(wb, ws_all, style_colhdr, rows = 2, cols = 1:ncols_all, gridExpand = TRUE)
setRowHeights(wb, ws_all, rows = 2, heights = 42)

# Data rows
alt <- FALSE
pval_cols_all <- which(grepl("^p \\(", colnames(sum_tbl)))
for (r in seq_len(nrow(sum_tbl))) {
  cur_row <- r + 2
  is_subhdr <- grepl("^  (Continuous|Binary)", sum_tbl[r, 1])
  writeData(wb, ws_all, as.data.frame(sum_tbl[r,]),
            startRow = cur_row, startCol = 1, colNames = FALSE)
  if (is_subhdr) {
    addStyle(wb, ws_all, style_subhdr, rows = cur_row,
             cols = 1:ncols_all, gridExpand = TRUE)
    mergeCells(wb, ws_all, cols = 1:ncols_all, rows = cur_row)
  } else {
    lbl_sty <- if (alt) style_altlbl else style_rowlbl
    dat_sty <- if (alt) style_alt    else style_data
    addStyle(wb, ws_all, lbl_sty, rows = cur_row, cols = 1, gridExpand = TRUE)
    addStyle(wb, ws_all, dat_sty, rows = cur_row, cols = 2:ncols_all, gridExpand = TRUE)

    # Highlight significant p-values in summary sheet
    for (pval_col in pval_cols_all) {
      pval_str <- as.character(sum_tbl[r, pval_col])
      if (!is.na(pval_str) && pval_str != "—" && pval_str != "") {
        pval_num <- suppressWarnings(as.numeric(gsub("<", "", pval_str)))
        if (!is.na(pval_num) && pval_num < 0.05) {
          addStyle(wb, ws_all, style_sig, rows = cur_row,
                   cols = pval_col, gridExpand = TRUE)
        }
      }
    }
    alt <- !alt
  }
  setRowHeights(wb, ws_all, rows = cur_row, heights = 16)
}
setColWidths(wb, ws_all, cols = 1, widths = 32)
setColWidths(wb, ws_all, cols = 2:ncols_all, widths = 13)
setColWidths(wb, ws_all, cols = pval_cols_all, widths = 10)
freezePane(wb, ws_all, firstActiveRow = 3, firstActiveCol = 2)

# ── Individual stratified sheets ───────────────────────────────────────────────
for (b_idx in seq_along(blocks)) {
  block  <- blocks[[b_idx]]
  ncols  <- ncol(block)
  sheet_name <- c("Race", "Obesity", "BreastCancer_F1", "Menarche")[b_idx]

  addWorksheet(wb, sheet_name)

  # Title row
  title_txt <- paste0("Table 1 — Intergen Sample Characteristics: ",
                       gsub("── Stratified by: ", "",
                            block[1, 1]))
  writeData(wb, sheet_name, data.frame(x = title_txt), startRow = 1,
            startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:ncols, rows = 1)
  addStyle(wb, sheet_name, style_title, rows = 1, cols = 1:ncols,
           gridExpand = TRUE)
  setRowHeights(wb, sheet_name, rows = 1, heights = 30)

  # Column headers
  writeData(wb, sheet_name,
            as.data.frame(t(colnames(block))),
            startRow = 2, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet_name, style_colhdr, rows = 2, cols = 1:ncols,
           gridExpand = TRUE)
  setRowHeights(wb, sheet_name, rows = 2, heights = 36)

  # Data rows
  data_start <- 3
  alt <- FALSE
  for (r in seq_len(nrow(block))) {
    cur_row <- data_start + r - 1
    cell_vals <- block[r, ]

    # Detect row type
    is_section <- grepl("^──", cell_vals[1, 1])
    is_subhdr  <- grepl("^  (Continuous|Binary)", cell_vals[1, 1])

    writeData(wb, sheet_name,
              as.data.frame(cell_vals),
              startRow = cur_row, startCol = 1, colNames = FALSE)

    if (is_section) {
      addStyle(wb, sheet_name, style_section, rows = cur_row,
               cols = 1:ncols, gridExpand = TRUE)
      mergeCells(wb, sheet_name, cols = 1:ncols, rows = cur_row)
      alt <- FALSE
    } else if (is_subhdr) {
      addStyle(wb, sheet_name, style_subhdr, rows = cur_row,
               cols = 1:ncols, gridExpand = TRUE)
      mergeCells(wb, sheet_name, cols = 1:ncols, rows = cur_row)
    } else {
      lbl_style  <- if (alt) style_altlbl  else style_rowlbl
      dat_style  <- if (alt) style_alt     else style_data

      addStyle(wb, sheet_name, lbl_style,  rows = cur_row,
               cols = 1,       gridExpand = TRUE)
      addStyle(wb, sheet_name, dat_style,  rows = cur_row,
               cols = 2:ncols, gridExpand = TRUE)

      # Highlight significant p-values in ALL p-value columns
      pval_cols <- which(grepl("^p[-( ]", colnames(block)))
      for (pval_col in pval_cols) {
        pval_str <- as.character(block[r, pval_col])
        if (!is.na(pval_str) && pval_str != "—" && pval_str != "") {
          pval_num <- suppressWarnings(as.numeric(gsub("<", "", pval_str)))
          if (!is.na(pval_num) && pval_num < 0.05) {
            addStyle(wb, sheet_name, style_sig, rows = cur_row,
                     cols = pval_col, gridExpand = TRUE)
          }
        }
      }
      alt <- !alt
    }
    setRowHeights(wb, sheet_name, rows = cur_row, heights = 16)
  }

  # Column widths — first col wide for labels, p-value cols narrower
  setColWidths(wb, sheet_name, cols = 1,       widths = 32)
  pval_col_idxs <- which(grepl("^p[-( ]", colnames(block)))
  data_col_idxs <- setdiff(2:ncols, pval_col_idxs)
  if (length(data_col_idxs) > 0)
    setColWidths(wb, sheet_name, cols = data_col_idxs, widths = 14)
  if (length(pval_col_idxs) > 0)
    setColWidths(wb, sheet_name, cols = pval_col_idxs, widths = 10)

  # Freeze panes
  freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 2)
}

# ── Footer note sheet ─────────────────────────────────────────────────────────
addWorksheet(wb, "Notes")
notes <- data.frame(Note = c(
  "Table 1 — Intergen WGBS Sample Characteristics",
  "",
  "Continuous variables: Mean (SD). Test: Wilcoxon rank-sum (2 groups) or Kruskal-Wallis (3+ groups).",
  "Binary variables: N (%). Test: Chi-squared with simulated p-value.",
  "Significant p-values (p < 0.05) are highlighted in dark red bold in the Excel version for review; remove highlighting if the target journal discourages color.",
  "",
  "Groups:",
  "  Race         — aakid: 0 = Non-Black, 1 = Black",
  "  Obesity      — HVobeseF1: 0 = Not Obese, 1 = Obese",
  "  Breast Ca F1 — breastcaF1: 0 = No, 1 = Yes",
  "  Menarche     — derived: youngmen_adol=1 → Young (<12), oldmen_adol=1 → Old (>14), else Middle",
  "",
  "ICE = Index of Concentration at the Extremes (Q4 = most advantaged quartile shown)"
))
writeData(wb, "Notes", notes, colNames = FALSE)
setColWidths(wb, "Notes", cols = 1, widths = 90)

# ── Save ──────────────────────────────────────────────────────────────────────
saveWorkbook(wb, opt$output, overwrite = TRUE)
message("Saved: ", opt$output)

# ── PDF output ────────────────────────────────────────────────────────────────
message("Generating PDF: ", opt$output_pdf)

# Helper to render one block as a grob (graphical table object)
render_block_pdf <- function(block, title_txt) {
  # Remove section/subheader rows that are just decorators
  is_section <- grepl("^──", block[, 1])
  is_subhdr  <- grepl("^  (Continuous|Binary|N$)", block[, 1])
  block_clean <- block[!(is_section | is_subhdr), , drop = FALSE]

  # Trim leading spaces from variable column for PDF
  block_clean[, 1] <- trimws(block_clean[, 1])

  # Build color vectors for rows
  n_rows  <- nrow(block_clean)
  fill_colors <- rep(c("white", "#EEF4FA"), length.out = n_rows)

  # Header fill
  header_fill <- rep("#2E75B6", ncol(block_clean))

  tt <- ttheme_minimal(
    core = list(
      bg_params = list(fill = fill_colors, col = "#CCCCCC", lwd = 0.5),
      fg_params = list(fontsize = 7, fontfamily = "sans", hjust = 0,
                       x = 0.04)
    ),
    colhead = list(
      bg_params = list(fill = "#2E75B6", col = "white", lwd = 0.5),
      fg_params = list(col = "white", fontsize = 7.5, fontface = "bold",
                       fontfamily = "sans")
    )
  )

  # Centre-align all columns except first
  n_cols <- ncol(block_clean)
  col_just <- c("left", rep("centre", n_cols - 1))

  tbl <- tableGrob(block_clean, rows = NULL, theme = tt)

  # Title grob
  title_grob <- textGrob(
    title_txt,
    gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "sans",
              col = "#1F4E79"),
    hjust = 0, x = 0.01
  )

  # Note grob
  note_grob <- textGrob(
    "Continuous: Mean (SD), Wilcoxon/Kruskal-Wallis. Binary: N (%), Chi-squared. *p<0.05 highlighted.",
    gp = gpar(fontsize = 6, fontfamily = "sans", col = "#555555", fontface = "italic"),
    hjust = 0, x = 0.01
  )

  arrangeGrob(
    title_grob,
    tbl,
    note_grob,
    ncol   = 1,
    heights = unit(c(0.6, nrow(block_clean) * 0.38 + 0.5, 0.4), "cm")
  )
}

block_titles <- c(
  "Table 1A — Stratified by Race (Black vs Non-Black)",
  "Table 1B — Stratified by Obesity (F1)",
  "Table 1C — Stratified by Breast Cancer (F1)",
  "Table 1D — Stratified by Menarche Timing"
)

pdf(opt$output_pdf, width = 22, height = 11, paper = "special")

for (b_idx in seq_along(blocks)) {
  block <- blocks[[b_idx]]

  # Remove purely decorative rows for PDF
  is_section <- grepl("^──", block[, 1])
  is_subhdr  <- grepl("^  (Continuous|Binary|N$)", block[, 1])
  block_clean <- block[!(is_section | is_subhdr), , drop = FALSE]
  block_clean[, 1] <- trimws(block_clean[, 1])

  # Mark significant p-values with asterisk in ALL p-value columns
  pval_cols <- which(grepl("^p[-( ]", colnames(block_clean)))
  for (pval_col in pval_cols) {
    for (r in seq_len(nrow(block_clean))) {
      pval_str <- as.character(block_clean[r, pval_col])
      if (!is.na(pval_str) && pval_str != "" && pval_str != "—") {
        pval_num <- suppressWarnings(as.numeric(gsub("<", "", pval_str)))
        if (!is.na(pval_num) && pval_num < 0.05) {
          block_clean[r, pval_col] <- paste0(pval_str, " *")
        }
      }
    }
  }

  n_rows <- nrow(block_clean)
  n_cols <- ncol(block_clean)

  # Row fill colors — alternating
  fill_mat <- matrix(
    rep(c("white", "#EEF4FA"), length.out = n_rows * n_cols),
    nrow = n_rows, ncol = n_cols, byrow = FALSE
  )
  # Highlight significant p-value cells red
  for (pval_col in pval_cols) {
    for (r in seq_len(n_rows)) {
      pval_str <- as.character(block_clean[r, pval_col])
      if (grepl("[*]", pval_str)) fill_mat[r, pval_col] <- "#FFE4E4"
    }
  }

  tt <- ttheme_minimal(
    core = list(
      bg_params  = list(fill = fill_mat, col = "#CCCCCC", lwd = 0.5),
      fg_params  = list(fontsize = 7, fontfamily = "sans"),
      padding    = unit(c(4, 3), "mm")
    ),
    colhead = list(
      bg_params  = list(fill = "#2E75B6", col = "white", lwd = 0.5),
      fg_params  = list(col = "white", fontsize = 7.5, fontface = "bold",
                        fontfamily = "sans"),
      padding    = unit(c(4, 3), "mm")
    )
  )

  tbl_grob <- tableGrob(block_clean, rows = NULL, theme = tt)

  # Set column widths: first col wide, p-value cols narrow, data cols medium
  n_cols_tbl    <- ncol(block_clean)
  pval_col_idxs <- which(grepl("^p[-( ]", colnames(block_clean)))
  col_w <- rep(2.8, n_cols_tbl)
  col_w[1] <- 5.0
  col_w[pval_col_idxs] <- 1.8
  tbl_grob$widths  <- unit(col_w, "cm")
  tbl_grob$heights <- unit(rep(0.52, nrow(tbl_grob)), "cm")

  title_grob <- textGrob(
    block_titles[b_idx],
    gp  = gpar(fontsize = 11, fontface = "bold", col = "#1F4E79"),
    hjust = 0, x = 0.01
  )

  note_grob <- textGrob(
    paste0("Continuous variables: Mean (SD) | Wilcoxon/Kruskal-Wallis test. ",
           "Binary variables: N (%) | Chi-squared test. ",
           "* p < 0.05"),
    gp  = gpar(fontsize = 6.5, col = "#555555", fontface = "italic"),
    hjust = 0, x = 0.01
  )

  # Calculate total table height including header row
  tbl_height_cm <- (n_rows + 1) * 0.52 + 0.5  # +1 for header row

  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))  # white background

  # Use viewport with margins
  pushViewport(viewport(
    x = 0.02, y = 0.02,
    width = 0.96, height = 0.96,
    just = c("left", "bottom")
  ))

  pushViewport(viewport(layout = grid.layout(
    nrow   = 3,
    heights = unit(c(0.9, tbl_height_cm, 0.5), "cm")
  )))

  # Title
  pushViewport(viewport(layout.pos.row = 1))
  grid.draw(title_grob)
  popViewport()

  # Table — draw inside its own viewport so it never overlaps title/note
  pushViewport(viewport(layout.pos.row = 2))
  pushViewport(viewport(x = 0, y = 1, just = c("left", "top"),
                        width = 1, height = 1))
  grid.draw(tbl_grob)
  popViewport()
  popViewport()

  # Note
  pushViewport(viewport(layout.pos.row = 3))
  grid.draw(note_grob)
  popViewport()

  popViewport()  # layout
  popViewport()  # margin viewport
}

dev.off()
message("Saved PDF: ", opt$output_pdf)