#!/usr/bin/env bash
# ============================================================
# run_12a_12b_hive.sh -- annotation (12a) + enrichment (12b) on Hive
#
# Steps 0-1 need internet  -> run on a LOGIN node (once).
# Steps 2-4 are offline     -> run on a compute node (sbatch / srun).
# Paths below are from your earlier run logs; adjust if they moved.
# ============================================================
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

PROJ=/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George
SCRIPTS=/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/new_not_tested           # folder holding 12a, 12b and helper.R (adjust)
REF=$PROJ/reference                                # caches live here
MOD_V2=$PROJ/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds
ANNOT_V2=$PROJ/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs
CPG=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/cpgIslandExt.txt.gz
mkdir -p "$REF"

STEP=${1:-help}
case "$STEP" in

# ---- 0) LOGIN NODE: check package versions + Enrichr reachability --------------
check)
  pixi run Rscript -e '
    for (p in c("rGREAT","annotatr","clusterProfiler","enrichR","msigdbr","org.Hs.eg.db","reactome.db","GO.db",
                "TxDb.Hsapiens.UCSC.hg38.knownGene","optparse","ChIPseeker"))
      cat(sprintf("%-36s %s\n", p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED"))
    cat("GREAT v4 hg38 table shipped with rGREAT:",
        nzchar(system.file("extdata","GREATv4.genes.hg38.tsv.gz", package="rGREAT")), "\n")
    cat("enrichR has background argument:", "background" %in% names(formals(enrichR::enrichr)), "\n")
    options(enrichR.base.address="https://maayanlab.cloud/Enrichr/", enrichR.live=TRUE, enrichR.quiet=TRUE)
    x <- tryCatch(enrichR::listEnrichrDbs(), error=function(e) e)
    cat("Enrichr reachable:", !inherits(x, "error") && !is.null(x), "\n")'
  ;;

# ---- 1) LOGIN NODE: build both caches (internet) --------------------------------
caches)
  pixi run Rscript $SCRIPTS/12a_annotate_modules.R --build_cache_only TRUE --genome hg38 \
    --annotatr_cache $REF/annotatr_hg38.rds
  pixi run Rscript $SCRIPTS/12b_enrichment.R --build_gene_set_cache_only TRUE \
    --gene_set_cache $REF/gene_sets_hg38.rds --gene_set_source native
  # If KEGG REST or MSigDB downloads fail, use one MSigDB release for all six instead:
  #   --gene_set_source msigdb
  ;;

# ---- 2) COMPUTE: 12a annotation (v2) -----------------------------------------------
annotate)
  pixi run Rscript $SCRIPTS/12a_annotate_modules.R \
    --project_root $PROJ --modules_v2 $MOD_V2 \
    --annotatr_cache $REF/annotatr_hg38.rds \
    --plot_modules yellow,brown,red,lightsteelblue1
  # defaults: --gene_assignment great --tss_source GREAT:hg38 --enrichment_background all
  ;;

# ---- 3) COMPUTE: 12b quick test on two modules ---------------------------------------
enrich_test)
  pixi run Rscript $SCRIPTS/12b_enrichment.R \
    --project_root $PROJ --annotation_dir_v1 $ANNOT_V2 \
    --enrichr_dbs "GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,MSigDB_Hallmark_2020,MSigDB_Oncogenic_Signatures,ClinVar_2025,WikiPathways_2024_Human,GTEx_Tissue_Expression_Up,Reactome_Pathways_2024,Human_Phenotype_Ontology,Epigenomics_Roadmap_HM_ChIP-seq,Elsevier_Pathway_Collection,Panther_2016,UK_Biobank_GWAS_v1,GWAS_Catalog_2025" \
    --gene_set_cache $REF/gene_sets_hg38.rds --modules yellow,brown,red,lightsteelblue1 \
    --plot_width 10 \
    --plot_height 7
  ;;

# ---- 4) COMPUTE (or a node with internet, for Enrichr): 12b, all modules -------------
enrich_all)
  pixi run Rscript $SCRIPTS/12b_enrichment.R \
    --project_root $PROJ --annotation_dir_v1 $ANNOT_V2 \
    --gene_set_cache $REF/gene_sets_hg38.rds
  # Enrichr needs internet; if this node has none, the log says
  # "ENRICHR WILL NOT RUN" and GREAT + clusterProfiler still complete.
  # Re-run Enrichr alone where there is internet with:
  #   --do_great FALSE --do_clusterprofiler FALSE
  ;;

*)
  echo "usage: $0 {check|caches|annotate|enrich_test|enrich_all}"
  ;;
esac
