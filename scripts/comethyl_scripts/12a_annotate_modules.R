#!/usr/bin/env Rscript

# ============================================================
# 12a_annotate_modules.R
#
# Purpose:
#   Annotate comethyl module regions with genic context, FANTOM5
#   enhancers and CpG context using annotatr -- the same annotation
#   codes comethyl::annotateModule() uses -- then test each module for
#   enrichment of those features and plot their distribution.
#
# What this script does:
#   - loads one or more module objects (v1/v2/v3 network variants)
#   - annotates ALL regions (grey included) with annotatr categories:
#       genic : promoter (<1 kb upstream of TSS), 1-5 kb upstream,
#               5' UTR, exon, intron, 3' UTR, intergenic
#       enhancer : FANTOM5 permissive enhancers (phase 1+2)
#       CpG   : island, shore (0-2 kb), shelf (2-4 kb), open sea
#     A region gets every label it overlaps (annotatr behaviour), plus
#     one "primary" label per family by precedence for stacked plots.
#   - assigns genes to regions the way comethyl::annotateModule() does:
#     GREAT basalPlusExt (5 kb up / 1 kb down of each TSS, extended up to
#     1 Mb), up to two genes per region ("A | B"), with distance to TSS,
#     run OFFLINE via rGREAT with GREAT v4's own hg38 gene table
#     (--tss_source GREAT:hg38, default) or TxDb knownGene
#     (--tss_source TxDb:hg38). Gene description / Ensembl / Entrez IDs
#     come from org.Hs.eg.db (comethyl used a live biomaRt query).
#     --gene_assignment nearest restores the old single nearest gene.
#   - Fisher's exact test per module x category against the chosen
#     background, BH-FDR within each family (genic+enhancer; CpG)
#   - optional width-adjusted logistic-regression sensitivity estimate
#   - writes annotated region tables, gene lists, background genes,
#     distribution tables, enrichment tables, plots and run logs
#
# What this script does NOT do:
#   - pathway/GO enrichment (Script 12b)
#   - module-trait significance filtering
#
# ------------------------------------------------------------
# annotatr downloads (CpG islands, FANTOM5 enhancers, liftOver chain)
# need internet. Hive compute nodes have none, so build the annotation
# cache ONCE on a LOGIN node:
#
#   pixi run Rscript 12a_annotate_modules.R \
#     --build_cache_only TRUE \
#     --genome hg38 \
#     --annotatr_cache /path/to/project/reference/annotatr_hg38.rds
#
# then point every compute-node run at that cache with --annotatr_cache.
# Fully offline alternative: also pass --cpg_island_file and
# --enhancer_file (+ --enhancer_chain_file) when building the cache.
#
# Typical run:
#   Rscript 12a_annotate_modules.R \
#     --project_root /path/to/project \
#     --modules_v2 /path/to/v2/Modules_for_downstream.rds \
#     --annotatr_cache /path/to/project/reference/annotatr_hg38.rds \
#     --plot_modules yellow,brown,red,lightsteelblue1
#
# Options (see --help for full text):
#   --project_root            project root; outputs go under comethyl_output/12_annotation/
#   --modules_v1/_v2/_v3      module RDS per network variant (at least one)
#   --sample_info             optional; only logged (kept for compatibility)
#   --sample_id_col           optional; only logged
#   --genome                  genome assembly (default hg38)
#   --helper_file             path to helper.R (default: next to this script)
#   --annotatr_cache          RDS of built annotatr annotations (load if exists, else build+save)
#   --rebuild_annotatr_cache  force rebuilding the cache (default FALSE)
#   --build_cache_only        build/save the cache and exit (run on a login node)
#   --cpg_island_file         local UCSC cpgIslandExt.txt.gz (offline CpG build)
#   --enhancer_file           local FANTOM5 enhancer BED (offline enhancer build)
#   --enhancer_chain_file     liftOver chain for --enhancer_file (e.g. hg19ToHg38.over.chain.gz)
#   --include_enhancers       annotate FANTOM5 enhancers (default TRUE)
#   --gene_assignment         great (default; comethyl's GREAT rule) | nearest (single nearest gene)
#   --tss_source              GREAT:hg38 (default; GREAT v4 genes) | TxDb:hg38 (knownGene protein-coding)
#   --great_mode              basalPlusExt (default) | twoClosest | oneClosest
#   --basal_upstream / --basal_downstream / --extension   GREAT domain sizes (5000 / 1000 / 1e6)
#   --enrichment_background   all (default; every region incl. grey) | assigned (non-grey only)
#   --adjust_width            add width-adjusted logistic-regression ORs (default FALSE; slower)
#   --plot_modules            comma list of modules for the "selected" plots + forest plots
#   --color_overrides         e.g. "yellow=darkgoldenrod3,white=grey80" for module bar colours
#   --log2OR_cap              colour-scale cap for heatmaps (default 3)
#   --fdr_cut                 FDR threshold for stars/filled points (default 0.05)
#
# Output structure per variant:
#   <project_root>/comethyl_output/12_annotation/<cpg_label>/<region_label>/<variant_name>/
#     Annotated_Regions.tsv / .xlsx   (gene_symbol / distance_to_TSS as "A | B")
#     Region_Gene_Associations.tsv    (one row per region-gene pair; used by 12b)
#     Module_Gene_List.tsv, Module_Gene_Summary.tsv / .xlsx
#     Background_Genes.txt / .tsv
#     Module_Annotation_Distribution.tsv / .xlsx
#     Module_Genic_Enrichment.tsv / .xlsx    (genic + enhancer categories)
#     Module_CpG_Enrichment.tsv / .xlsx
#     plots/  (PDF via cairo_pdf + 600 dpi PNG)
#       Genic_Enrichment_Heatmap, CpG_Enrichment_Heatmap
#       Genic_Distribution_{stacked,overlap}_all, CpG_Distribution_{stacked,overlap}_all
#       *_selected versions and Genic/CpG_Enrichment_Forest_selected (if --plot_modules)
#     run_parameters.txt, run_log.txt, session_info.txt
# ============================================================

message("Starting Script 12a")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(openxlsx)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(AnnotationDbi)
  library(stringr)
  library(ggplot2)
})

# ============================================================
# 1) CLI
# ============================================================
option_list <- list(
  make_option("--project_root", type = "character", default = NULL,
              help = "Project root directory. Outputs are written under <project_root>/comethyl_output/12_annotation/. Required unless --build_cache_only TRUE."),
  make_option("--modules_v1", type = "character", default = NULL,
              help = "Module RDS for variant v1 (all PCs). Must contain $regions with RegionID, chr, start, end, module."),
  make_option("--modules_v2", type = "character", default = NULL,
              help = "Module RDS for variant v2 (exclude protected PCs)."),
  make_option("--modules_v3", type = "character", default = NULL,
              help = "Module RDS for variant v3 (technical PCs only)."),
  make_option("--sample_info", type = "character", default = NULL,
              help = "Optional sample info file (xlsx/csv/tsv). Not used by annotation; only its dimensions and checksum are logged."),
  make_option("--sample_id_col", type = "character", default = NULL,
              help = "Optional sample ID column name in --sample_info (logged only)."),
  make_option("--genome", type = "character", default = "hg38",
              help = "Genome assembly with annotatr support [default %default]."),
  make_option("--helper_file", type = "character", default = NULL,
              help = "Path to helper.R [default: helper.R next to this script, else working directory]."),
  make_option("--annotatr_cache", type = "character", default = NULL,
              help = "RDS of built annotatr annotations. Loaded if it exists; otherwise annotations are built and saved here. Strongly recommended on HPC."),
  make_option("--rebuild_annotatr_cache", type = "character", default = "FALSE",
              help = "TRUE to rebuild --annotatr_cache even if it exists [default %default]."),
  make_option("--build_cache_only", type = "character", default = "FALSE",
              help = "TRUE to build/save --annotatr_cache and exit without annotating modules. Run on a node with internet unless local files are given [default %default]."),
  make_option("--cpg_island_file", type = "character", default = NULL,
              help = "Local UCSC cpgIslandExt.txt(.gz) for the genome. If given, CpG features are built from it with annotatr's algorithm instead of downloading."),
  make_option("--enhancer_file", type = "character", default = NULL,
              help = "Local FANTOM5 enhancer BED (e.g. human_permissive_enhancers_phase_1_and_2.bed.gz, hg19). If given, used instead of downloading."),
  make_option("--enhancer_chain_file", type = "character", default = NULL,
              help = "liftOver chain applied to --enhancer_file (e.g. hg19ToHg38.over.chain.gz). Omit if the BED is already in --genome."),
  make_option("--include_enhancers", type = "character", default = "TRUE",
              help = "TRUE to annotate and test FANTOM5 enhancers [default %default]."),
  make_option("--gene_assignment", type = "character", default = "great",
              help = "How genes are assigned to regions: 'great' = comethyl's GREAT rule (up to two genes per region, offline via rGREAT); 'nearest' = single nearest gene body [default %default]."),
  make_option("--tss_source", type = "character", default = "GREAT:hg38",
              help = "Gene/TSS set for --gene_assignment great: 'GREAT:hg38' = GREAT v4's own gene table shipped with rGREAT (matches comethyl); 'TxDb:hg38' = TxDb.Hsapiens.UCSC.hg38.knownGene protein-coding genes [default %default]."),
  make_option("--great_mode", type = "character", default = "basalPlusExt",
              help = "GREAT association rule: basalPlusExt (comethyl default), twoClosest or oneClosest [default %default]."),
  make_option("--basal_upstream", type = "double", default = 5000,
              help = "basalPlusExt basal domain upstream of TSS, bp [default %default]."),
  make_option("--basal_downstream", type = "double", default = 1000,
              help = "basalPlusExt basal domain downstream of TSS, bp [default %default]."),
  make_option("--extension", type = "double", default = 1e6,
              help = "Maximum domain extension, bp [default %default]."),
  make_option("--enrichment_background", type = "character", default = "all",
              help = "Background for Fisher's tests: 'all' = every region passed into network construction, grey included; 'assigned' = non-grey regions only (module vs rest of network) [default %default]."),
  make_option("--adjust_width", type = "character", default = "FALSE",
              help = "TRUE to add width-adjusted ORs from logistic regression (overlap ~ in_module + log(width)) as a sensitivity analysis [default %default]."),
  make_option("--plot_modules", type = "character", default = NULL,
              help = "Comma-separated modules of interest, e.g. 'yellow,brown,red,lightsteelblue1'. Produces *_selected distribution plots and forest plots."),
  make_option("--color_overrides", type = "character", default = NULL,
              help = "Module bar colour overrides, e.g. 'yellow=darkgoldenrod3,white=grey80'."),
  make_option("--log2OR_cap", type = "double", default = 3,
              help = "Cap for the heatmap log2(OR) colour scale; values beyond are squished [default %default]."),
  make_option("--fdr_cut", type = "double", default = 0.05,
              help = "FDR threshold for significance marks in plots [default %default].")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Annotate comethyl modules with annotatr genic/enhancer/CpG context and test enrichment."))

to_bool <- function(x, nm) {
  v <- toupper(trimws(as.character(x)))
  if (v %in% c("TRUE", "T", "1", "YES")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO")) return(FALSE)
  stop("--", nm, " must be TRUE or FALSE (got '", x, "')", call. = FALSE)
}
opt$rebuild_annotatr_cache <- to_bool(opt$rebuild_annotatr_cache, "rebuild_annotatr_cache")
opt$build_cache_only       <- to_bool(opt$build_cache_only, "build_cache_only")
opt$include_enhancers      <- to_bool(opt$include_enhancers, "include_enhancers")
opt$adjust_width           <- to_bool(opt$adjust_width, "adjust_width")

# ============================================================
# 2) Helper
# ============================================================
script_dir <- tryCatch({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
helper_file <- if (!is.null(opt$helper_file)) opt$helper_file else file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) {
  stop("helper.R not found at: ", helper_file, "\nPass --helper_file /path/to/helper.R.", call. = FALSE)
}
source(helper_file)
if (exists("setup_publication_font", mode = "function")) setup_publication_font()

# ============================================================
# 3) Validation
# ============================================================
opt$enrichment_background <- match.arg(opt$enrichment_background, c("all", "assigned"))
opt$gene_assignment <- match.arg(opt$gene_assignment, c("great", "nearest"))
opt$great_mode <- match.arg(opt$great_mode, c("basalPlusExt", "twoClosest", "oneClosest"))
if (opt$gene_assignment == "great") invisible(great_parse_tss_source(opt$tss_source, opt$genome))
for (nm in c("cpg_island_file", "enhancer_file", "enhancer_chain_file")) {
  if (!is.null(opt[[nm]])) validate_file_exists(opt[[nm]], paste0("--", nm))
}
if (!is.null(opt$enhancer_chain_file) && is.null(opt$enhancer_file)) {
  stop("--enhancer_chain_file given without --enhancer_file.", call. = FALSE)
}
if (opt$build_cache_only && is.null(opt$annotatr_cache)) {
  stop("--build_cache_only TRUE requires --annotatr_cache.", call. = FALSE)
}
if (!is.finite(opt$log2OR_cap) || opt$log2OR_cap <= 0) stop("--log2OR_cap must be > 0.", call. = FALSE)
if (!is.finite(opt$fdr_cut) || opt$fdr_cut <= 0 || opt$fdr_cut >= 1) stop("--fdr_cut must be in (0, 1).", call. = FALSE)

plot_modules <- if (!is.null(opt$plot_modules)) trimws(strsplit(opt$plot_modules, ",")[[1]]) else NULL
color_overrides <- annot_parse_color_overrides(opt$color_overrides)

if (!opt$build_cache_only) {
  stop_if_missing(opt$project_root, "--project_root")
  if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root, call. = FALSE)
  if (all(vapply(opt[c("modules_v1", "modules_v2", "modules_v3")], is.null, TRUE))) {
    stop("Provide at least one of --modules_v1/--modules_v2/--modules_v3.", call. = FALSE)
  }
  for (nm in c("modules_v1", "modules_v2", "modules_v3")) {
    if (!is.null(opt[[nm]])) validate_file_exists(opt[[nm]], paste0("--", nm))
  }
  if (!is.null(opt$sample_info)) validate_file_exists(opt$sample_info, "--sample_info")
  AnnotationHub::setAnnotationHubOption("CACHE", value = file.path(opt$project_root, ".cache"))
}

# ============================================================
# 4) Annotations (built once, shared across variants)
# ============================================================
annotations <- build_annotatr_annotations(
  genome = opt$genome,
  include_enhancers = opt$include_enhancers,
  cache_file = opt$annotatr_cache,
  rebuild = opt$rebuild_annotatr_cache,
  cpg_island_file = opt$cpg_island_file,
  enhancer_file = opt$enhancer_file,
  enhancer_chain_file = opt$enhancer_chain_file,
  verbose = TRUE
)
annot_prov <- attr(annotations, "provenance")
if (is.null(annot_prov)) annot_prov <- list(note = "no provenance attribute on cached annotations")

if (opt$build_cache_only) {
  message("Annotation cache ready: ", opt$annotatr_cache)
  message(paste(names(annot_prov), unlist(annot_prov), sep = ": ", collapse = "\n"))
  message("Script 12a: cache build finished")
  quit(save = "no", status = 0)
}

# GREAT regulatory domains (built once, shared across variants)
extended_tss <- NULL
if (opt$gene_assignment == "great") {
  extended_tss <- great_build_extended_tss(opt$tss_source, genome = opt$genome, mode = opt$great_mode,
                                           basal_upstream = opt$basal_upstream,
                                           basal_downstream = opt$basal_downstream,
                                           extension = opt$extension, verbose = TRUE)
}

sample_info_df <- NULL
if (!is.null(opt$sample_info)) {
  ext <- tolower(tools::file_ext(opt$sample_info))
  sample_info_df <- switch(ext,
    xlsx = openxlsx::read.xlsx(opt$sample_info), xls = openxlsx::read.xlsx(opt$sample_info),
    csv = utils::read.csv(opt$sample_info, check.names = FALSE),
    utils::read.delim(opt$sample_info, check.names = FALSE))
}

# ============================================================
# 5) Per-variant runner
# ============================================================
derive_out_dir <- function(modules_rds, project_root, step_name = "12_annotation") {
  variant_dir <- dirname(modules_rds); region_dir <- dirname(variant_dir)
  info <- list(variant_name = basename(variant_dir), region_label = basename(region_dir),
               cpg_label = basename(dirname(region_dir)),
               pipeline_root = file.path(project_root, "comethyl_output"))
  info$step_dir <- file.path(info$pipeline_root, step_name)
  info$out_dir <- file.path(info$step_dir, info$cpg_label, info$region_label, info$variant_name)
  safe_dir_create(info$out_dir)
  info
}

write_tsv_xlsx <- function(df, stem, xlsx = TRUE) {
  utils::write.table(df, paste0(stem, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (xlsx) openxlsx::write.xlsx(df, paste0(stem, ".xlsx"), rowNames = FALSE)
}

run_annotation_for_variant <- function(modules_path, variant_label) {
  di <- derive_out_dir(modules_path, opt$project_root)
  out_dir <- di$out_dir
  plot_dir <- file.path(out_dir, "plots"); safe_dir_create(plot_dir)
  log_file <- file.path(out_dir, "run_log.txt")
  L <- function(...) append_log(log_file, ...)

  L("Starting annotation for variant: ", variant_label, " (annotatr engine)")
  L("modules_path: ", modules_path)
  L("out_dir: ", out_dir)
  L("genome: ", opt$genome, "; enrichment_background: ", opt$enrichment_background,
    "; include_enhancers: ", opt$include_enhancers, "; adjust_width: ", opt$adjust_width)
  for (nm in names(annot_prov)) L("annotations ", nm, ": ", annot_prov[[nm]])
  if (!is.null(sample_info_df)) L("sample_info rows: ", nrow(sample_info_df), "; cols: ", ncol(sample_info_df))

  regions <- extract_regions_from_module_object(readRDS(modules_path), label = modules_path)
  L("Loaded regions: ", nrow(regions), "; modules (incl. grey): ", length(unique(regions$module)))

  # --- annotatr genic / enhancer / CpG context ---
  annotated <- annotate_regions_annotatr(regions, annotations, genome = opt$genome, verbose = TRUE)

  # --- gene assignment ---
  if (opt$gene_assignment == "great") {
    assoc <- great_region_gene_associations(regions, extended_tss)
    ginfo <- tryCatch(gene_info_from_orgdb(unique(assoc$gene_entrezID), opt$genome),
                      error = function(e) { L("WARNING: gene descriptions not added: ", conditionMessage(e)); NULL })
    collapsed <- great_collapse_associations(assoc, regions, ginfo)
    annotated <- cbind(annotated, collapsed[, setdiff(colnames(collapsed), "RegionID")])
    gp <- attr(extended_tss, "great_params")
    L("Gene assignment: GREAT ", gp$mode, " (", gp$tss_source, "; basal ", gp$basal_upstream, " up / ",
      gp$basal_downstream, " down; extension ", gp$extension, "; rGREAT ", gp$rGREAT_version, "; ",
      gp$n_genes, " genes)")
    tab <- table(factor(pmin(collapsed$n_genes, 3), levels = 0:3, labels = c("0", "1", "2", "3+")))
    L("Genes per region: ", paste(names(tab), tab, sep = "=", collapse = "; "))
  } else {
    gr <- GenomicRanges::GRanges(regions$chr, IRanges::IRanges(regions$start, regions$end))
    ng <- offline_nearest_gene(gr, genome = opt$genome, verbose = TRUE)
    ng <- ng[!duplicated(ng[, c("chr", "start", "end")]), ]
    annotated <- dplyr::left_join(annotated, ng[, c("chr", "start", "end", "gene_symbol", "gene_entrezID")],
                                  by = c("chr", "start", "end"))
    assoc <- annotated[!is.na(annotated$gene_symbol) & annotated$gene_symbol != "",
                       c("RegionID", "chr", "start", "end", "gene_entrezID", "gene_symbol", "module")]
    assoc$distance_to_TSS <- NA_integer_
    L("Gene assignment: single nearest gene (offline_nearest_gene)")
  }
  stopifnot(nrow(annotated) == nrow(regions))
  L("Regions with gene_context: ", sum(!is.na(annotated$gene_context)),
    "; with CpG_context: ", sum(!is.na(annotated$CpG_context)),
    "; with >=1 gene: ", length(unique(assoc$RegionID)))

  write_tsv_xlsx(annotated, file.path(out_dir, "Annotated_Regions"))
  utils::write.table(assoc[, c("RegionID", "chr", "start", "end", "module", "gene_symbol",
                               "gene_entrezID", "distance_to_TSS")],
                     file.path(out_dir, "Region_Gene_Associations.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # --- gene lists / background genes (one row per module-gene; "A | B" split) ---
  gl <- assoc %>% dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::distinct(module, gene_symbol) %>% dplyr::arrange(module, gene_symbol)
  utils::write.table(gl, file.path(out_dir, "Module_Gene_List.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  gs <- annotated %>% dplyr::count(module, name = "n_regions") %>%
    dplyr::left_join(gl %>% dplyr::count(module, name = "n_unique_genes"), by = "module") %>%
    dplyr::mutate(n_unique_genes = dplyr::coalesce(n_unique_genes, 0L)) %>%
    dplyr::arrange(dplyr::desc(n_regions), module)
  write_tsv_xlsx(gs, file.path(out_dir, "Module_Gene_Summary"))
  bg_genes <- sort(unique(gl$gene_symbol))
  write_lines_safe(bg_genes, file.path(out_dir, "Background_Genes.txt"))
  utils::write.table(data.frame(gene_symbol = bg_genes), file.path(out_dir, "Background_Genes.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  # --- distribution + enrichment ---
  dist <- summarizeAnnotationDistribution(annotated)
  write_tsv_xlsx(dist, file.path(out_dir, "Module_Annotation_Distribution"))

  enr <- list()
  for (fam in c("genic", "cpg")) {
    enr[[fam]] <- moduleAnnotatrEnrichment(annotated, family = fam,
                                           background = opt$enrichment_background,
                                           adjust_width = opt$adjust_width, verbose = TRUE)
    stem <- file.path(out_dir, if (fam == "genic") "Module_Genic_Enrichment" else "Module_CpG_Enrichment")
    write_tsv_xlsx(enr[[fam]], stem)
    L("Wrote ", basename(stem), " (", nrow(enr[[fam]]), " tests; ",
      sum(enr[[fam]]$fdr < opt$fdr_cut), " at FDR < ", opt$fdr_cut, ")")
  }

  # --- plots ---
  fam_name <- c(genic = "Genic", cpg = "CpG")
  for (fam in c("genic", "cpg")) {
    plotAnnotatrEnrichmentHeatmap(enr[[fam]], p_cut = opt$fdr_cut, log2OR_cap = opt$log2OR_cap,
                                  title = paste0(fam_name[[fam]], " enrichment by module (", variant_label, ")"),
                                  file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Enrichment_Heatmap")))
    if (opt$adjust_width) {
      plotAnnotatrEnrichmentHeatmap(enr[[fam]], p_cut = opt$fdr_cut, log2OR_cap = opt$log2OR_cap,
                                    use_width_adjusted = TRUE,
                                    file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Enrichment_Heatmap_width_adjusted")))
    }
    for (mode in c("stacked", "overlap")) {
      plotAnnotationDistribution(dist, family = fam, mode = mode, color_overrides = color_overrides,
                                 file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Distribution_", mode, "_all")))
      if (!is.null(plot_modules)) {
        plotAnnotationDistribution(dist, family = fam, mode = mode, groups = plot_modules,
                                   color_overrides = color_overrides,
                                   file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Distribution_", mode, "_selected")))
      }
    }
    if (!is.null(plot_modules)) {
      plotAnnotatrEnrichmentForest(enr[[fam]], groups = plot_modules, p_cut = opt$fdr_cut,
                                   color_overrides = color_overrides,
                                   file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Enrichment_Forest_selected")))
    }
  }
  L("Wrote plots to: ", plot_dir)

  # --- provenance ---
  params <- c(
    paste0("timestamp\t", timestamp_now()),
    paste0("annotation_engine\tannotatr"),
    paste0("variant_label\t", variant_label),
    paste0(names(opt), "\t", vapply(opt, function(x) if (is.null(x)) "" else paste(x, collapse = ","), "")),
    paste0("annotations_", names(annot_prov), "\t", unlist(annot_prov)),
    paste0("n_regions\t", nrow(regions)),
    paste0("n_modules_incl_grey\t", length(unique(regions$module))),
    paste0("n_background_genes\t", length(bg_genes)),
    if (!is.null(extended_tss)) paste0("great_", names(attr(extended_tss, "great_params")), "\t",
                                       unlist(attr(extended_tss, "great_params"))) else character(0)
  )
  write_lines_safe(params, file.path(out_dir, "run_parameters.txt"))
  write_session_info(file.path(out_dir, "session_info.txt"),
                     extra_files = list(modules_path = modules_path, sample_info = opt$sample_info,
                                        annotatr_cache = opt$annotatr_cache,
                                        cpg_island_file = opt$cpg_island_file,
                                        enhancer_file = opt$enhancer_file))
  L("Finished annotation for variant: ", variant_label)
  invisible(out_dir)
}

# ============================================================
# 6) Run variants
# ============================================================
variants <- c(modules_v1 = "v1_all_pcs", modules_v2 = "v2_exclude_protected_pcs",
              modules_v3 = "v3_technical_pcs_only")
for (nm in names(variants)) {
  if (!is.null(opt[[nm]])) run_annotation_for_variant(opt[[nm]], variants[[nm]])
}

message("Script 12a complete: Annotate Modules finished")
