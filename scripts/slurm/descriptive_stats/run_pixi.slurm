#!/bin/bash
#SBATCH --job-name=Cometh_Intergen_CHDS
#SBATCH --mail-user=gekuodza@health.ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=23:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=100G
#SBATCH --account=lasallegrp #genome-center-grp, publicgrp ,lasallegrp
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

# # Sanity check
# pixi run Rscript --vanilla -e 'cat("R version:", R.version.string, "\n"); library(comethyl); cat("comethyl OK\n")'

##############
# Run Script #
##############
# First run the metadata cleaning scripts to ensure that the metadata is in the correct format for downstream analysis. 
# Then, run the table generation scripts to create summary tables and trait correlation analyses.
pixi run /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/clean_metadata/00_rename_metadata_publication_labels.R

pixi run /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/clean_metadata/01_check_intergen_priority_variables.R

#Table generation scripts
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/descriptive_stats/01_table1.R \
  --input /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/numeric_conversion_check/Intergen173_070702026_full_dataset_priority_numeric.xlsx \
  --out_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/descriptive_stats/01_table1 \
  --traits   /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/table1_traits.txt \
  --na_codes "F,N"

# pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/tables/01_table1.R \
#   --input /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_publication_labels.xlsx \
#   --out_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/00_table1 \
#   --traits   /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/table1_traits.txt \
#   --na_codes "F,N"

#Correlation analysis for traits

pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/descriptive_stats/02_trait_correlations.R \
  --input   /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_publication_labels.xlsx \
  --traits  /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/00_figure1_heatmap_primary_traits.txt \
  --out_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/descriptive_stats/02_trait_correlations \
  --na_codes "F,N"

# pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/descriptive_stats/02_trait_correlations.R \
#   --input    /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/numeric_conversion_check/Intergen173_070702026_full_dataset_priority_numeric.xlsx  \
#   --traits   /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/table1_traits.txt \
#   --out_dir  /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/descriptive_stats/02_trait_correlations \
#   --na_codes "F,N"

#random script to rename metadata for publication
# pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/clean_metadata/rename_metadata_publication_labels.R \
#   --input_path /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026.xlsx \
#   --output_path /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_publication_labels.xlsx \
#   --dictionary_path "/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Data Dictionary for Intergen Class File (3).xlsx" \
#   --dictionary_output_path /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_DataDictionary_publication_labels.xlsx

end=$(date +%s)
echo "Run completed in $((end - start)) seconds"
date