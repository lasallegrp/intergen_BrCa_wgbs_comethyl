#!/bin/bash
#SBATCH --job-name=Cometh_Intergen_CHDS
#SBATCH --mail-user=gekuodza@health.ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=23:00:00
#SBATCH --partition=low
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=250G
#SBATCH --account=publicgrp #genome-center-grp, publicgrp ,lasallegrp
#SBATCH --output=00_std_err_logs/comethyl_pearson_report.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report.err
#SBATCH --chdir=/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George

set -euo pipefail

###################
# Run Information #
###################
start=$(date +%s)
hostname
date

####################
# Load Environment #
####################
source ~/.profile

# Pixi settings for cluster stability
export PIXI_HOME=/quobyte/lasallegrp/George/.pixi
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME

mkdir -p "$PIXI_CACHE_DIR"

# Keep R isolated from user libraries/configs
unset R_LIBS_USER
export R_LIBS_USER=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# Optional threading hygiene
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Move into Pixi workspace
cd /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/reproducibility/pixi

conda activate Comethyl_v1.3.0

# # Sanity check
# Rscript --vanilla -e 'cat("R version:", R.version.string, "\n"); library(comethyl); cat("comethyl OK\n")'

##############
# Run Script #
##############

# ##############
# # Comethyl Scripts#
# ##############
# # # Script 00
# # Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/00_import_cpg_reports.R \
# #   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
# #   --meta_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_publication_labels.xlsx \
# #   --cyto_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/processed/cytosine_reports \
# #   --workers 4

# # # Script 01
# # Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/01_filter_cpgs_call_regions.R \
# #   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
# #   --input_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/00_cpg_extraction/Unfiltered_BSseq.rds \
# #   --cov 3 \
# #   --per_sample 0.75

# # #Script 02
# # Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/02_filter_regions.R \
# #   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
# #   --input_regions /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/01_filter_regions/cov3_75pct/Regions.txt \
# #   --cov_min 4 \
# #   --meth_sd 0.08


# # # Script 03
# # Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/03_region_methylation.R \
# #   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
# #   --input_bs /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/01_filter_regions/cov3_75pct/Filtered_BSseq.rds \
# #   --input_regions /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/02_filter_regions/cov3_75pct/covMin10_methSD0p07/Filtered_Regions.txt


# Script 04
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/04_pcs_and_trait_diagnostics.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
#   --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_full_dataset_priority_numeric.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
#   --threads 4


# # Script 04
# #v1
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/04_pcs_and_trait_diagnostics.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
#   --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_full_dataset_priority_numeric.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
#   --pc_model intercept_only \
#   --primary_cor_method pearson \
#   --threads 4

# #v2
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/04_pcs_and_trait_diagnostics.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_full_dataset_priority_numeric.xlsx \
  --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
  --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
  --pc_model protected \
  --primary_cor_method bicor \
  --threads 4


# #old version
# # Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/04_pcs_and_trait_diagnostics.R \
# #   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
# #   --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
# #   --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/CHDS_CHARGE_combined_010926.xlsx \
# #   --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
# #   --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
# #   --primary_cor_method bicor \
# #   --threads 4

# # Script 05 bicor
# # Script 05A: v1 original comethyl-style adjustment
# # PC derivation = intercept_only
# # Adjustment = all selected PCs

# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/05_adjust_region_methylation.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
#   --input_pcs /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/intercept_only/PCs.rds \
#   --pc_trait_stats /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/intercept_only/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
#   --run_v1 TRUE \
#   --run_v2 FALSE \
#   --run_v3 FALSE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05 \
#   --v3_p_thresh 0.05

# # Script 05B: v2 protected-PC exclusion using bicor
# # PC derivation = protected_model
# # Adjustment = exclude protected-associated PCs, then adjust remaining PCs

# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/05_adjust_region_methylation.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
#   --input_pcs /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/protected_model/PCs.rds \
#   --pc_trait_stats /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/protected_model/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
#   --run_v1 FALSE \
#   --run_v2 TRUE \
#   --run_v3 FALSE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05 \
#   --v3_p_thresh 0.05



# # Script 05
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/05_adjust_region_methylation.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
#   --input_pcs /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/PCs.rds \
#   --pc_trait_stats /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
#   --run_v1 TRUE \
#   --run_v2 TRUE \
#   --run_v3 FALSE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05

# # Script 06
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/06_soft_power.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --input_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
#   --input_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --threads 4


# # Script 07
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/07_module_detection.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --regions_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/02_filter_regions/cov3_75pct/covMin10_methSD0p07/Filtered_Regions.txt \
#   --input_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
#   --softpower_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/06_soft_power/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/SoftPower_pearson.rds \
#   --input_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --softpower_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/06_soft_power/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/SoftPower_pearson.rds \
#   --module_cor pearson \
#   --deep_split 4 \
#   --min_module_size 35 \
#   --merge_cut_height 0.1 \
#   --threads 4

# # Script 08
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/08_module__sample_diagnostics.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --modules_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/07_module_detection/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules.rds \
#   --modules_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/07_module_detection/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules.rds \
#   --module_cor pearson \
#   --sample_dendro_distance euclidean \
#   --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_publication_labels.xlsx


#script 09a
Rscript --vanilla /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/09a_me_trait_analysis.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_full_dataset_priority_numeric.xlsx \
  --modules_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
  --modules_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --run_bicor TRUE \
  --p_thresh 0.05 \
  --top_n 25 \
  --module_dendro_distance bicor 

#script 09b
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/09b_me_trait_presentation_heatmaps.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --stats_file_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --modules_rds_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --set_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets \
  --p_thresh 0.05 \
  --top_n 25 \
  --top_width   14 \
  --top_height  7   \
  --min_sig_traits 1 \
  --focused_width 15 \
  --focused_height 11 \
  --module_order_mode dendro \
  --module_order_trait_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/module_ordering/exposure_outcome_traits.tsv \
  --module_order_p_thresh 0.05 \
  --module_order_score p \
  --trait_order_mode input \
  --drop_grey_modules TRUE


  # --module_order_mode exposure_outcome \
    # --trait_order_mode exposure_outcome \
# # #Script 10
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/10_circos_and_module_summaries.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --modules_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
#   --modules_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
#   --genome hg38 \
#   --single_module_color red \
#   --all_modules_width 4.0 \
#   --all_modules_height 4.0 \
#   --single_module_width 5.5 \
#   --single_module_height 5.5 \
#   --counts_width 7 \
#   --counts_height 6


# #Script 11
#v1
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/11_me_trait_pair_plots.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --modules_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
  --stats_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_full_dataset_priority_numeric.xlsx \
  --set_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets \
  --p_thresh 0.05 \
  --only_significant TRUE \
  --cat_plot_style violin \
  --top_n_pairs 50

#v2
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/11_me_trait_pair_plots.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --modules_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --stats_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_full_dataset_priority_numeric.xlsx\
  --set_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets \
  --p_thresh 0.05 \
  --only_significant TRUE \
  --cat_plot_style violin \
  --continuous_width 8 \
  --continuous_height 7 \
  --categorical_width 8 \
  --categorical_height 7 \
  --top_n_pairs 50

#script 12a
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/12a_annotate_modules.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_thesis.xlsx \
  --modules_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
  --modules_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --genome hg38 \
  --annotation_mode offline \
  --helper_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/helper.R \
  --cpg_island_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/cpgIslandExt.txt.gz


# #script 12b
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/12b_enrichment.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --annotation_dir_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v1_all_pcs \
#   --annotation_dir_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs \
#   --trait_list_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets \
#   --stats_file_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --stats_file_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --p_thresh 0.05 \
#   --do_kegg false \
#   --do_go false \
#   --do_reactome false \
#   --do_enrichr true \
#   --enrichr_dbs GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,MSigDB_Hallmark_2020,MSigDB_Oncogenic_Signatures,ClinVar_2025,WikiPathways_2024_Human \
#   --use_enrichr_background false \
#   --enrichr_retries 5 \
#   --enrichr_sleep 2 \
#   --plot_p_col pvalue \
#   --plot_top_n 10 \
#   --plot_sig_cutoff 0.05

#  #full list
#   #  --enrichr_dbs ClinVar_2025,WikiPathways_2024_Human,DisGeNET,GTEx_Tissue_Expression_Up,Human_Phenotype_Ontology \
# #  --trait_list_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/12_heatmap_primary_traits.txt \

#script 12b just one text file only
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/12b_enrichment.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --annotation_dir_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v1_all_pcs \
  --annotation_dir_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs \
  --stats_file_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --stats_file_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --module_list_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/enrichrment_modules/modules.txt \
  --p_thresh 0.05 \
  --do_kegg true \
  --do_go true \
  --do_reactome false \
  --do_enrichr true \
  --enrichr_dbs GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,MSigDB_Hallmark_2020,MSigDB_Oncogenic_Signatures,ClinVar_2025,WikiPathways_2024_Human,GTEx_Tissue_Expression_Up,Reactome_Pathways_2024,Human_Phenotype_Ontology,Epigenomics_Roadmap_HM_ChIP-seq,Elsevier_Pathway_Collection,Panther_2016,UK_Biobank_GWAS_v1,GWAS_Catalog_2025 \
  --use_enrichr_background false \
  --enrichr_retries 5 \
  --enrichr_sleep 2 \
  --plot_p_col p.adjust \
  --plot_top_n 10 \
  --plot_sig_cutoff 0.05 \
  --plot_width 10 \
  --plot_height 7


# #script 12b just one text file only(Sylvia thesis)
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/12b_enrichment.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --annotation_dir_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v1_all_pcs \
#   --annotation_dir_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs \
#   --stats_file_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --stats_file_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --module_list_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/enrichrment_modules/modules.txt \
#   --p_thresh 0.05 \
#   --do_kegg false \
#   --do_go false \
#   --do_reactome false \
#   --do_enrichr true \
#   --enrichr_dbs GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,GWAS_Catalog_2025 \
#   --use_enrichr_background false \
#   --enrichr_retries 5 \
#   --enrichr_sleep 2 \
#   --plot_p_col pvalue \
#   --plot_top_n 8 \
#   --plot_sig_cutoff 0.05 \
#   --plot_width 10 \
#   --plot_height 8

# #p.adjust", "pvalue

# # Script 13
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/13_run_module_cancer_gene_enrichment.R \
#   --cancer_file "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Census_allMon May 25 17_26_07 2026.csv" \
#   --annot_file "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Annotated_Regions.xlsx" \
#   --output_step_name "13_cancer_gene_enrichment" \
#   --module_col "module" \
#   --gene_col "gene_symbol" \
#   --module_trait_file "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv" \
#   --trait_module_col "module" \
#   --trait_col "trait" \
#   --trait_p_col "p" \
#   --trait_effect_col "bicor" \
#   --trait_fdr_cutoff 1 \
#   --cancer_keyword "breast" \
#   --top_cancer_types 20 \
#   --plot_width 12 \
#   --plot_height 8 \
#   --dpi 600

# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/13_run_module_cancer_gene_enrichment.R \
#   --cancer_file "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Census_allMon May 25 17_26_07 2026.csv" \
#   --annot_file "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Annotated_Regions.xlsx" \
#   --out_dir "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/13_cancer_gene_enrichment/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs" \
#   --output_step_name "13_cancer_gene_enrichment" \
#   --module_col "module" \
#   --gene_col "gene_symbol" \
#   --cancer_fdr_cutoff 0.05 \
#   --module_trait_file "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/11_module_trait/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Module_Trait_Associations.csv" \
#   --trait_module_col "Module" \
#   --trait_col "trait" \
#   --trait_fdr_col "NONE" \
#   --trait_p_col "p" \
#   --trait_effect_col "cor" \
#   --trait_fdr_cutoff 0.05 \
#   --cancer_keyword "breast" \
#   --top_cancer_types 20 \
#   --focus_modules "yellow,brown,tan,salmon,magenta,darkgrey" \
#   --plot_width 12 \
#   --plot_height 8 \
#   --dpi 600

# # Script 13a
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/13A_core_regulatory_state_overlap.R \
#   --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
#   --region_version covMin10_methSD0p07 \
#   --annotation_dir_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/covMin10_methSD0p07/v2_exclude_protected_pcs \
#   --modules turquoise \
#   --source roadmap 

  # --annotation_dir_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/covMin10_methSD0p07/v1_all_pcs \

# # Script 14
# Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/tables/00_table1.R \
#   --input /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173.xlsx \
#   --output /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/00_table1/table1_intergen.xlsx \
#   --output_pdf /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/00_table1/table1_intergen.pdf \
#   --verify

#random script to rename metadata for publication
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/clean_metadata/rename_metadata_publication_labels.R \
  --input_path /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026.xlsx \
  --output_path /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_publication_labels.xlsx \
  --dictionary_path "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Data Dictionary for Intergen Class File (3).xlsx" \
  --dictionary_output_path /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_DataDictionary_publication_labels.xlsx

end=$(date +%s)
echo "Run completed in $((end - start)) seconds"
date