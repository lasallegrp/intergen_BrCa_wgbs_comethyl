# =============================================================================
# check_PCs_vs_celltype.R
#
# Question: are the methylation PCs I protected for race (PC1, PC2) also
#           capturing blood cell-type composition?
#
# What this script does
#   Step 1  Load the Comethyl PCs, IDOL cell proportions, and sample traits
#   Step 2  Match all three tables by sample ID
#   Step 3  Correlate each PC with each cell type (a table you can read)
#   Step 4  Does neutrophil fraction differ by race?
#   Step 5  For PC1 and PC2: how much is explained by race, by cells, by both?
#
# Run:  pixi run Rscript check_PCs_vs_celltype.R
# =============================================================================

library(data.table)   # fread() / fwrite() for reading and writing tables

# ==== SETTINGS YOU CAN MODIFY ================================================
pcs_file         <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/PCs.rds"                    # MODIFY: your saved Comethyl PCs
cell_prop_file   <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/cell_composition_chr22/cellprop_idol_CP_normalized.tsv"
trait_file       <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_Class07_07_2026.xlsx"          # MODIFY: table with sample ID + race
trait_id_column  <- "LabID"                              # MODIFY: name of the sample-ID column
race_column      <- "aakid"                          # MODIFY: 1 = Black, 0 = not (or a factor)
n_pcs_to_check   <- 33                                     # look at PC1..PC5, not only 1 and 2
output_directory <- "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/pc_celltype_check"
# =============================================================================

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)


# ---- Step 1: Load the three inputs -----------------------------------------
# Input : PCs.rds        -> matrix, rows = samples, columns = PC1, PC2, ...
#         cell proportions-> one row per sample: sample, CD8T, CD4T, NK, Bcell, Mono, Neu
#         traits          -> one row per sample, includes the race column
# Output: three R objects with the same information

# readRDS() loads a single R object saved earlier with saveRDS().
pc_matrix <- readRDS(pcs_file)

# colClasses = list(character = "sample") keeps IDs like "4001" as text,
# so they match the row names of pc_matrix (which are always text).
cell_proportions <- fread(cell_prop_file, colClasses = list(character = "sample"))
# read_excel() reads an .xlsx sheet (readxl package).
   #   sheet = 1  -> first sheet; MODIFY to the sheet name, e.g. sheet = "traits"
   # as.data.table() converts the result so the rest of the script works unchanged.
   sample_traits <- as.data.table(readxl::read_excel(trait_file, sheet = 1))

   # Excel stores LabID as a number (4042); turn it into text "4042" so it
   # matches the sample IDs in the PC and cell-proportion tables.
   sample_traits[[trait_id_column]] <- as.character(sample_traits[[trait_id_column]])

# Quick look -- healthy result: ~173 rows each, sample IDs look like "4001"
print(dim(pc_matrix)); print(head(rownames(pc_matrix)))
print(dim(cell_proportions)); print(dim(sample_traits))


# ---- Step 2: Match samples across the three tables -------------------------
# Input : the three objects above, possibly in different sample orders
# Output: one table, analysis_table, one row per sample with PCs + cells + race

pc_table <- data.table(sample = rownames(pc_matrix),
                       pc_matrix[, 1:n_pcs_to_check])

# merge() joins two tables on a shared column; by = tells it which one.
# Only samples present in BOTH tables are kept.
analysis_table <- merge(pc_table, cell_proportions, by = "sample")
analysis_table <- merge(analysis_table,
                        sample_traits[, c(trait_id_column, race_column), with = FALSE],
                        by.x = "sample", by.y = trait_id_column)

# Neutrophil-to-lymphocyte ratio on the log scale (a single "myeloid vs
# lymphoid" number; log makes it roughly symmetric).
analysis_table$log_NLR <- log(analysis_table$Neu /
                              (analysis_table$CD4T + analysis_table$CD8T +
                               analysis_table$NK + analysis_table$Bcell))

cat("Samples matched in all three tables:", nrow(analysis_table), "\n")
# Healthy: 173. If much lower, the sample IDs are formatted differently.


# ---- Step 3: Correlate each PC with each cell type --------------------------
# Input : analysis_table
# Output: pc_cell_correlations.tsv -- rows = PCs, columns = cell types (Pearson r)
#
# How to read it: |r| > 0.5 means that PC is largely a cell-composition axis.
colnames(readRDS(pcs_file))[1:5]
pc_names   <- paste0("PC_", 1:n_pcs_to_check)         # "PC_1", "PC_2", ...
cell_names <- c("Neu", "CD4T", "CD8T", "Bcell", "NK", "Mono", "log_NLR")

# cor(x, y) with two tables gives every column of x against every column of y.
pc_cell_correlations <- cor(analysis_table[, ..pc_names],
                            analysis_table[, ..cell_names])
# Note: "..pc_names" tells data.table "use the column names stored in pc_names".

print(round(pc_cell_correlations, 2))
fwrite(data.table(PC = pc_names, round(pc_cell_correlations, 3)),
       file.path(output_directory, "pc_cell_correlations.tsv"), sep = "\t")


# ---- Step 4: Does neutrophil fraction differ by race? -----------------------
# Input : analysis_table
# Output: printed group means and a t-test
#
# Why: if race and neutrophil fraction are related, then "variation associated
# with race" and "variation associated with cells" overlap, and protecting one
# partly protects the other.

race_values <- analysis_table[[race_column]]   # [[ ]] pulls one column out by name
print(tapply(analysis_table$Neu, race_values, mean))   # mean Neu per race group
print(t.test(Neu ~ get(race_column), data = analysis_table))
# get(race_column) turns the text "race_black" into the actual column.


# ---- Step 5: For each protected PC, separate race from cell composition ----
# Input : analysis_table
# Output: pc_race_vs_cells.tsv with, for each PC:
#           r2_race        - variance explained by race alone
#           r2_cells       - variance explained by Neu + log_NLR alone
#           r2_both        - by race and cells together
#           race_p_adjusted- is race still associated after accounting for cells?
#
# lm(y ~ x) fits a linear regression; summary(fit)$r.squared is the fraction
# of the PC's variation the predictors explain (0 to 1).

results_list <- list()
for (pc in pc_names[1:2]) {                              # MODIFY: which PCs you protected
  formula_race  <- as.formula(paste(pc, "~", race_column))
  formula_cells <- as.formula(paste(pc, "~ Neu + log_NLR"))
  formula_both  <- as.formula(paste(pc, "~", race_column, "+ Neu + log_NLR"))

  fit_race  <- lm(formula_race,  data = analysis_table)
  fit_cells <- lm(formula_cells, data = analysis_table)
  fit_both  <- lm(formula_both,  data = analysis_table)

  # coef(summary(fit)) is the coefficient table; row 2 is the race term,
  # column 4 is its p-value.
  results_list[[pc]] <- data.table(
    PC              = pc,
    r2_race         = summary(fit_race)$r.squared,
    r2_cells        = summary(fit_cells)$r.squared,
    r2_both         = summary(fit_both)$r.squared,
    race_p_alone    = coef(summary(fit_race))[2, 4],
    race_p_adjusted = coef(summary(fit_both))[2, 4]
  )
}
pc_race_vs_cells <- rbindlist(results_list)   # stack the per-PC rows into one table
print(pc_race_vs_cells)
fwrite(pc_race_vs_cells, file.path(output_directory, "pc_race_vs_cells.tsv"), sep = "\t")

# How to read pc_race_vs_cells:
#   r2_cells small (< 0.1)                    -> PC is not a cell-composition axis. Good.
#   r2_cells large, race_p_adjusted still small -> PC carries BOTH race and cell signal.
#   r2_cells large, race_p_adjusted large       -> race association runs through cell
#                                                 composition; protecting this PC mostly
#                                                 protects cell-type variation.