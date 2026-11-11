#!/usr/bin/env Rscript

# =============================================================================
# Rename Intergen metadata to publication-readable variable labels
#
# This script intentionally uses spaces, punctuation, and units in column names.
# These labels are suitable for figures and tables. In downstream R code, refer
# to them with backticks, e.g. metadata$`DDT `.
#
# Usage:
# Rscript rename_metadata_publication_labels.R \
#   --input_path /path/Intergen173_070702026_thesis.xlsx \
#   --output_path /path/Intergen173_070702026_publication_labels.xlsx \
#   --dictionary_path "/path/Data Dictionary for Intergen Class File.xlsx" \
#   --dictionary_output_path /path/Intergen173_070702026_DataDictionary_publication_labels.xlsx

  #  META=/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata

  #  pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/clean_metadata/00_rename_metadata_publication_labels.R \
  #    --input_path  $META/Intergen173_full_dataset_priority_numeric_with_IDOL_cellprop.xlsx \
  #    --output_path $META/Intergen173_full_dataset_priority_numeric_with_IDOL_publication_labels.xlsx \
  #    --dictionary_output_path $META/Intergen173_DataDictionary_with_IDOL_publication_labels.xlsx
# ==============================================================================

suppressPackageStartupMessages(library(openxlsx))

parse_args <- function(args) {
  if (length(args) == 0L) return(list())
  if (length(args) %% 2L != 0L) {
    stop("Arguments must be supplied as --name value pairs.")
  }
  keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]
  stats::setNames(as.list(vals), keys)
}

`%||%` <- function(x, y) if (is.null(x) || !nzchar(x)) y else x

args <- parse_args(commandArgs(trailingOnly = TRUE))

input_path <- args$input_path %||%
  "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_Class07_07_2026.xlsx"

output_path <- args$output_path %||%
  file.path(dirname(input_path), "Intergen173_070702026_publication_labels.xlsx")

dictionary_path <- args$dictionary_path %||%
  "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Data_Dictionary_for_Intergen_Multiomics_Class_File_07072026.xlsx"

dictionary_output_path <- args$dictionary_output_path %||%
  file.path(dirname(output_path), "Intergen173_070702026_DataDictionary_publication_labels.xlsx")

if (!file.exists(input_path)) stop("Input metadata file not found: ", input_path)
if (!file.exists(dictionary_path)) stop("Input data dictionary not found: ", dictionary_path)

# -----------------------------------------------------------------------------
# Publication-readable labels: raw metadata variable -> display-ready label
# -----------------------------------------------------------------------------
rename_map <- c(
  # Technical / batch variables
  "batch2" = "Batch",
  "order2" = "Order within batch",
  "batchorder2" = "Batch order",

  # Participant characteristics
  "F1sex" = "F1 sex",
  "aakid" = "F1 Black race",
  "black_race" = "F1 Black race",
  "Cati_Age" = "F1 age (years)",
  "HVobeseF1" = "F1 obesity",
  "HVoverwtF1" = "F1 overweight",

  # Perinatal DDT / DDE
  "ddt" = "DDT",
  "ddtq1" = "DDT: Q1",
  "ddtq2" = "DDT: Q2",
  "ddtq3" = "DDT: Q3",
  "ddtq4" = "DDT: Q4",
  "ddtT1" = "DDT: T1",
  "ddtT2" = "DDT: T2",
  "ddtT3" = "DDT: T3",
  "dde" = "DDE",
  "ddeq1" = "DDE: Q1",
  "ddeq2" = "DDE: Q2",
  "ddeq3" = "DDE: Q3",
  "ddeq4" = "DDE: Q4",
  "ddeT1" = "DDE: T1",
  "ddeT2" = "DDE: T2",
  "ddeT3" = "DDE: T3",
  "op_ddt" = "op_DDT",

  # ICE indices
  "ICE60edu1_q1" = "ICE education: Q1",
  "ICE60edu1_q2" = "ICE education: Q2",
  "ICE60edu1_q3" = "ICE education: Q3",
  "ICE60edu1_q4" = "ICE education: Q4",
  "ICE60rac1_q1" = "ICE race: Q1",
  "ICE60rac1_q2" = "ICE race: Q2",
  "ICE60rac1_q3" = "ICE race: Q3",
  "ICE60rac1_q4" = "ICE race: Q4",
  "ICE60inc1_q1" = "ICE income: Q1",
  "ICE60inc1_q2" = "ICE income: Q2",
  "ICE60inc1_q3" = "ICE income: Q3",
  "ICE60inc1_q4" = "ICE income: Q4",
  "ICE60incxrac1_q1" = "ICE income-race: Q1",
  "ICE60incxrac1_q2" = "ICE income-race: Q2",
  "ICE60incxrac1_q3" = "ICE income-race: Q3",
  "ICE60incxrac1_q4" = "ICE income-race: Q4",

  # Reproductive / breast cancer risk factors
  "agemenarche_adol" = "F1 age at menarche",
  "youngmen_adol" = "F1 early menarche",
  "oldmen_adol" = "F1 late menarche",
  "breastcaF1" = "F1 Breast cancer",
  "breastcaF0" = "F0 Breast cancer",

  # Birth characteristics
  "ga_wks" = "F1 gestational age (weeks)",
  "under37" = "F1 preterm birth",
  "iugrz" = "F1 IUGR",

  # Menopausal status at midlife visit
  "postmenopausal" = "F1 postmenopausal",
  "perimenopausal" = "F1 perimenopausal",

  # Perinatal organochlorines / PFAS metabolites
  "EtPFOSAAcOH" = "EtPFOSAAcOH",
  "HCB" = "HCB",
  "MePFOSAAcOH" = "MePFOSAAcOH",
  "Oxy" = "Oxychlordane",

  # Perinatal PCBs
  "PCB28" = "PCB28",
  "PCB66" = "PCB66",
  "PCB74" = "PCB74",
  "PCB99" = "PCB99",
  "PCB101" = "PCB101",
  "PCB105" = "PCB105",
  "PCB118" = "PCB118",
  "PCB138" = "PCB138",
  "PCB153" = "PCB153",
  "PCB156" = "PCB156",
  "PCB167" = "PCB167",
  "PCB170" = "PCB170",
  "PCB180" = "PCB180",
  "PCB183" = "PCB183",
  "PCB187" = "PCB187",
  "PCB194" = "PCB194",
  "PCB201" = "PCB201",
  "PCB203" = "PCB203",
  "pcb15" = "PCB15",
  "pcb56" = "PCB56",
  "pcb82" = "PCB82",
  "pcb146" = "PCB146",
  "pcb174" = "PCB174",
  "pcb177" = "PCB177",
  "pcb178" = "PCB178",
  "pcb199" = "PCB199",

  # Perinatal PFAS
  "PFBS" = "PFBS",
  "PFDeA" = "PFDeA",
  "PFDoA" = "PFDoA",
  "PFHpA" = "PFHpA",
  "PFHxA" = "PFHxA",
  "PFHxS" = "PFHxS",
  "PFNA" = "PFNA",
  "PFOA" = "PFOA",
  "PFOS" = "PFOS",
  "PFOSA" = "PFOSA",
  "PFUdA" = "PFUdA",

  # Perinatal lipids
  "TC" = "F1 Total cholesterol",
  "TG" = "F1 Triglycerides",
  "TL" = "F1 Total lipids",
  "bBHC" = "F1 bBHC",

  # Neurodegeneration biomarkers
  #the variables measure by ALZ path are the same as those measured by MSD, but the ALZ path assay is more sensitive and specific for p-tau217. 
  #The ALZ path assay is used in the main analyses, but the MSD assay is included here for completeness.
  "GFAP__pg_mL_" = "F1 GFAP",
  "NFL__pg_mL_" = "F1 NfL",
  "Ab40__pg_mL_" = "F1 AB40",
  "Ab42__pg_mL_" = "F1 AB42",
  "Ab42_Ab40_ratio" = "F1 AB42/AB40 ratio",
  "Total_Tau__pg_mL_" = "F1 Total tau",
  "MSD_pTau_217__pg_mL_" = "F1 p-tau217, MSD",
  "MSD_pTau_217__total_tau_ratio" = "F1 p-tau217/total tau, MSD",
  "MSD_pTau_217__A_42_ratio" = "F1 p-tau217/AB42, MSD",
  "Total_Tau__A_42_ratio" = "F1 Total tau/AB42 ratio",
  "ALZ_pTau_217__pg_mL_" = "F1 p-tau217",
  "ALZ_pTau_217__total_tau_ratio" = "F1 p-tau217/total tau",
  "ALZ_pTau_217__A_42_ratio" = "F1 p-tau217/AB42",

  # Kidney function
  "Creatinine__cr__mg_dL_" = "F1 Creatinine (mg/dL)",
  "Cystatin_C____cys__mg_L_" = "F1 Cystatin C (mg/L)",

  # APOE
  "APOE_Genotype" = "F1 APOE genotype",
  "APOe4ext" = "F1 APOE e4 carrier",

  # Allostatic load
  "AL_SUM10" = "F1 Allostatic load score",
  "AL_CRP" = "F1 Allostatic load: high CRP",
  "AL_DHEAS" = "F1 Allostatic load: low DHEA-S",
  "AL_cholratio" = "F1 Allostatic load: high total/HDL ratio",
  "AL_HDL" = "F1 Allostatic load: low HDL",
  "AL_IL6" = "F1 Allostatic load: high IL-6",
  "AL_HBA1C" = "F1 Allostatic load: high HbA1c",
  "AL_BPDI" = "F1 Allostatic load: high DBP",
  "AL_BPSYS" = "F1 Allostatic load: high SBP",
  "AL_PBF" = "F1 Allostatic load: high body fat",
  "AL_WC" = "F1 Allostatic load: high waist circumference",

  # IDOL cell type proportions
  "IDOL_Neu"     = "Neutrophils",
  "IDOL_CD4T"    = "CD4+ T cells",
  "IDOL_CD8T"    = "CD8+ T cells",
  "IDOL_Bcell"   = "B cells",
  "IDOL_NK"      = "NK cells",
  "IDOL_Mono"    = "Monocytes",
  "IDOL_log_NLR" = "log NLR",

  # Midlife cardiometabolic biomarkers
  "hsCRP" = "F1 hsCRP",
  "DHEAS" = "F1 DHEA-S",
  "cholratio" = "F1 Total/HDL cholesterol ratio",
  "HDLC" = "F1 HDL cholesterol",
  "IL6" = "F1 IL-6",
  "HbA1c" = "F1 HbA1c (%)",
  "bpdias_avg" = "F1 Diastolic blood pressure",
  "bpsys_avg" = "F1 Systolic blood pressure",
  "bodyfat_avg" = "F1 Body fat (%)",
  "bmi" = "F1 BMI",
  "waistcirc_avg" = "F1 Waist circumference",
  "weight_avg" = "F1 Weight (kg)"
)

# -----------------------------------------------------------------------------
# Read metadata. The first Excel column is the participant identifier / row name.
# -----------------------------------------------------------------------------
message("Reading metadata: ", input_path)
metadata <- read.xlsx(input_path, rowNames = TRUE, check.names = FALSE, sep.names = " ")
original_names <- colnames(metadata)

missing_from_map <- setdiff(original_names, names(rename_map))
if (length(missing_from_map) > 0L) {
  warning(
    "These variables are not in rename_map and will remain unchanged:\n  ",
    paste(missing_from_map, collapse = "\n  ")
  )
}

new_names <- unname(ifelse(
  original_names %in% names(rename_map),
  rename_map[original_names],
  original_names
))

if (anyDuplicated(new_names)) {
  duplicate_labels <- unique(new_names[duplicated(new_names)])
  stop(
    "Duplicate publication labels would be created:\n  ",
    paste(duplicate_labels, collapse = "\n  "),
    "\nResolve these before writing the data."
  )
}

colnames(metadata) <- new_names
message("Renamed ", sum(original_names != new_names), " of ", length(original_names), " variables.")

# -----------------------------------------------------------------------------
# Write renamed metadata without replacing the source dataset.
# -----------------------------------------------------------------------------
wb_data <- createWorkbook()
addWorksheet(wb_data, "data")
writeData(wb_data, "data", cbind(LabID = rownames(metadata), metadata), keepNA = FALSE)
freezePane(wb_data, "data", firstRow = TRUE, firstCol = TRUE)
setColWidths(wb_data, "data", cols = 1:(ncol(metadata) + 1), widths = "auto")
saveWorkbook(wb_data, output_path, overwrite = TRUE)
message("Wrote publication-label metadata: ", output_path)

# -----------------------------------------------------------------------------
# Build a revised data dictionary using the supplied source data dictionary.
# The source dictionary is expected to have Variable in column 2 and Label in 3.
# -----------------------------------------------------------------------------
source_dictionary <- read.xlsx(dictionary_path, sheet = 1, startRow = 2,
                               colNames = TRUE, check.names = FALSE)

if (!all(c("Variable", "Label") %in% colnames(source_dictionary))) {
  stop("Could not find columns named 'Variable' and 'Label' in: ", dictionary_path)
}

source_dictionary <- source_dictionary[!is.na(source_dictionary$Variable), , drop = FALSE]
source_label_map <- stats::setNames(as.character(source_dictionary$Label),
                                    as.character(source_dictionary$Variable))

revised_dictionary <- data.frame(
  Number = seq_along(original_names),
  Original_variable = original_names,
  Publication_label = new_names,
  Original_description = unname(source_label_map[original_names]),
  Renamed = original_names != new_names,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

unmatched_descriptions <- is.na(revised_dictionary$Original_description)
if (any(unmatched_descriptions)) {
  revised_dictionary$Original_description[unmatched_descriptions] <- ""
}

wb_dict <- createWorkbook()
addWorksheet(wb_dict, "dictionary")

header_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#2E4057", halign = "left",
  textDecoration = "bold", border = "Bottom"
)
body_style <- createStyle(halign = "left", valign = "top", wrapText = TRUE)
alternate_style <- createStyle(fgFill = "#F0F4F8")

writeData(wb_dict, "dictionary", revised_dictionary, headerStyle = header_style)
addStyle(wb_dict, "dictionary", body_style,
         rows = 2:(nrow(revised_dictionary) + 1),
         cols = 1:ncol(revised_dictionary), gridExpand = TRUE)
addStyle(wb_dict, "dictionary", alternate_style,
         rows = seq(2, nrow(revised_dictionary) + 1, by = 2),
         cols = 1:ncol(revised_dictionary), gridExpand = TRUE, stack = TRUE)
freezePane(wb_dict, "dictionary", firstRow = TRUE)
setColWidths(wb_dict, "dictionary", cols = 1, widths = 8)
setColWidths(wb_dict, "dictionary", cols = 2, widths = 35)
setColWidths(wb_dict, "dictionary", cols = 3, widths = 40)
setColWidths(wb_dict, "dictionary", cols = 4, widths = 90)
setColWidths(wb_dict, "dictionary", cols = 5, widths = 10)

saveWorkbook(wb_dict, dictionary_output_path, overwrite = TRUE)
message("Wrote revised data dictionary: ", dictionary_output_path)

# Save variables that need a future map entry, if any.
if (length(missing_from_map) > 0L) {
  unmapped_path <- file.path(dirname(output_path), "unmapped_metadata_variables.csv")
  write.csv(data.frame(variable = missing_from_map), unmapped_path, row.names = FALSE)
  message("Wrote unmapped-variable report: ", unmapped_path)
}

message("Done.")
