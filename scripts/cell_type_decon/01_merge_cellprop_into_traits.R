# =============================================================================
# merge_cellprop_into_traits.R
#
# Goal: add IDOL blood cell-type proportions (and log NLR) to the Intergen173
#       trait/metadata spreadsheet, so the Comethyl adjustment step can use them.
#
# Steps
#   Step 1  Read the trait spreadsheet (.xlsx) and the cell-proportion table (.tsv)
#   Step 2  Make the sample IDs comparable (both as text, e.g. "4042")
#   Step 3  Prefix cell columns with "IDOL_" and add log_NLR
#   Step 4  Merge by sample ID and check every sample matched
#   Step 5  Save a NEW file (the original spreadsheet is never overwritten)
#
# Run:  pixi run Rscript merge_cellprop_into_traits.R
# Needs: pixi add r-readxl r-writexl
# =============================================================================

library(data.table)   # fread(), merge(), fwrite()

# ==== SETTINGS YOU CAN MODIFY ================================================
project_dir <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George"

trait_file     <- file.path(project_dir, "data/metadata/numeric_conversion_check",
                            "Intergen173_070702026_full_dataset_priority_numeric.xlsx")
trait_sheet    <- 1            # MODIFY: sheet number or name, e.g. "Sheet1"
trait_id_column <- "rowname"   # the sample-ID column in the spreadsheet

cell_prop_file <- file.path(project_dir, "results/cell_composition_chr22",
                            "cellprop_idol_CP_normalized.tsv")
cell_id_column <- "sample"     # the sample-ID column in the cell-proportion table

output_dir  <- file.path(project_dir, "data/metadata")
output_name <- "Intergen173_full_dataset_priority_numeric_with_IDOL_cellprop"  # MODIFY
# =============================================================================


# ---- Step 1: Read both files ------------------------------------------------
# Input : trait .xlsx  -> one row per sample, ~137 columns (rowname, Batch, ...)
#         cell .tsv    -> one row per sample: sample, CD8T, CD4T, NK, Bcell, Mono, Neu
# Output: sample_traits, cell_proportions (data.tables)

# read_excel() reads one sheet of an .xlsx file (readxl package).
# as.data.table() converts it so we can use data.table tools below.
sample_traits <- as.data.table(readxl::read_excel(trait_file, sheet = trait_sheet))

# colClasses keeps IDs as text, so "4042" is not turned into the number 4042.
cell_proportions <- fread(cell_prop_file,
                          colClasses = list(character = cell_id_column))

cat("Trait table:", nrow(sample_traits), "rows,", ncol(sample_traits), "columns\n")
cat("Cell table :", nrow(cell_proportions), "rows,", ncol(cell_proportions), "columns\n")
# Healthy: 173 rows in both.


# ---- Step 2: Make sample IDs comparable --------------------------------------
# Input : the two ID columns (Excel may store "rowname" as a number)
# Output: both ID columns as plain text, e.g. "4042"
#
# as.integer() first drops any decimal ("4042.0" -> 4042), then as.character()
# turns it into text. trimws() removes stray spaces around the ID.
sample_traits[[trait_id_column]] <- trimws(as.character(as.integer(sample_traits[[trait_id_column]])))
cell_proportions[[cell_id_column]] <- trimws(cell_proportions[[cell_id_column]])

print(head(sample_traits[[trait_id_column]]))     # expect "4042" "4084" ...
print(head(cell_proportions[[cell_id_column]]))   # expect "4001" "4002" ...


# ---- Step 3: Rename cell columns and add log NLR ----------------------------
# Input : cell_proportions with columns CD8T, CD4T, NK, Bcell, Mono, Neu
# Output: same table with columns IDOL_CD8T ... IDOL_Neu, IDOL_log_NLR
#
# Why the prefix: makes it obvious where these columns came from, and avoids
# clashing with any existing column of the same name in the trait file.

cell_type_columns <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")

# log NLR = log( neutrophils / all lymphocytes ). One number for the
# myeloid-vs-lymphoid axis that drove PC1.
lymphocyte_total <- cell_proportions$CD8T + cell_proportions$CD4T +
                    cell_proportions$NK   + cell_proportions$Bcell
cell_proportions$log_NLR <- log(cell_proportions$Neu / lymphocyte_total)

# setnames(table, old, new) renames columns in place (data.table).
setnames(cell_proportions,
         old = c(cell_type_columns, "log_NLR"),
         new = paste0("IDOL_", c(cell_type_columns, "log_NLR")))

print(head(cell_proportions, 3))


# ---- Step 4: Merge and check ------------------------------------------------
# Input : sample_traits, cell_proportions
# Output: traits_with_cells -- the trait table plus 7 new IDOL_ columns
#
# merge(..., all.x = TRUE) keeps EVERY row of the trait table, even if a sample
# has no cell estimate (it gets NA). That way nothing silently disappears.

ids_only_in_traits <- setdiff(sample_traits[[trait_id_column]], cell_proportions[[cell_id_column]])
ids_only_in_cells  <- setdiff(cell_proportions[[cell_id_column]], sample_traits[[trait_id_column]])
cat("IDs in traits but not cell table:", length(ids_only_in_traits), "\n")
cat("IDs in cell table but not traits:", length(ids_only_in_cells), "\n")
if (length(ids_only_in_traits)) print(ids_only_in_traits)
if (length(ids_only_in_cells))  print(ids_only_in_cells)
# Healthy: 0 and 0.

traits_with_cells <- merge(sample_traits, cell_proportions,
                           by.x = trait_id_column, by.y = cell_id_column,
                           all.x = TRUE, sort = FALSE)   # sort = FALSE keeps the original row order

# stopifnot() halts with an error if the condition is FALSE -- a safety check.
stopifnot(nrow(traits_with_cells) == nrow(sample_traits))
cat("Samples with missing IDOL values:", sum(is.na(traits_with_cells$IDOL_Neu)), "\n")
# Healthy: 0.


# ---- Step 5: Save as a new file ---------------------------------------------
# Input : traits_with_cells
# Output: <output_name>.xlsx (for your Comethyl scripts that read Excel)
#         <output_name>.tsv  (plain text, easy to inspect / diff / version)

xlsx_path <- file.path(output_dir, paste0(output_name, ".xlsx"))
tsv_path  <- file.path(output_dir, paste0(output_name, ".tsv"))

writexl::write_xlsx(traits_with_cells, xlsx_path)   # writes one sheet
fwrite(traits_with_cells, tsv_path, sep = "\t")

cat("Wrote:\n ", xlsx_path, "\n ", tsv_path, "\n")
cat("Final table:", nrow(traits_with_cells), "rows,", ncol(traits_with_cells), "columns\n")
# Healthy: 173 rows, original column count + 7.