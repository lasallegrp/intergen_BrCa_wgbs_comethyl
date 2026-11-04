#!/usr/bin/env Rscript

# ============================================================
# 12b_enrichment.R
#
# Purpose:
#   Functional enrichment of comethyl modules, three ways:
#
#   1. GREAT (primary; region-based)  -- rGREAT local mode, offline.
#      Foreground = a module's regions; background = all regions passed
#      into network construction (grey included) by default. Same gene/TSS
#      set and association rule as Script 12a (GREAT v4 genes,
#      basalPlusExt), so the genes behind each term match 12a's gene
#      columns. Binomial test over region midpoints (expected = fraction
#      of background covered by the term's regulatory domains) plus
#      rGREAT's gene-level hypergeometric test.
#   2. clusterProfiler ORA (gene-based) -- enricher() on the same cached
#      gene sets, gene list = 12a's GREAT-assigned genes per module,
#      universe = 12a's Background_Genes (genes assigned to any region).
#   3. Enrichr (gene-based, online)     -- same gene lists, background =
#      Background_Genes (Speedrichr endpoint; needs enrichR >= 3.3 and
#      internet). One connectivity check at start; failures are logged
#      loudly and recorded in enrichr_status.tsv.
#
# Gene-set collections (GREAT + clusterProfiler):
#   GO_BP, GO_MF, GO_CC, KEGG, Reactome, Hallmark   (default: all)
#   Read from a gene-set cache RDS built ONCE on a node with internet:
#     Rscript 12b_enrichment.R --build_gene_set_cache_only TRUE \
#       --gene_set_cache /path/reference/gene_sets_hg38.rds [--gene_set_source native|msigdb]
#   native (default): GO from org.Hs.eg.db, KEGG from KEGG REST,
#                     Reactome from reactome.db, Hallmark from msigdbr
#   msigdb:           all six from one MSigDB release via msigdbr
#
# Module selection (default: ALL non-grey modules):
#   --modules yellow,brown      explicit comma list, or
#   --module_list_file          .txt (or directory of .txt; one run per file), or
#   --trait_list_file + --stats_file_v*   modules with trait p <= --p_thresh
#
# Options (see --help for full text):
#   --project_root, --annotation_dir_v1 (required; 12a output folder), --annotation_dir_v2/_v3
#   --modules, --module_list_file, --trait_list_file, --stats_file_v1/_v2/_v3, --p_thresh
#   --collections               comma list or 'all' [all]
#   --gene_set_cache            gene-set cache RDS (required unless building it)
#   --build_gene_set_cache_only / --gene_set_source
#   --do_great / --do_clusterprofiler / --do_enrichr   [TRUE / TRUE / TRUE]
#   --tss_source, --great_mode, --basal_upstream, --basal_downstream, --extension
#                               GREAT settings (default: comethyl's; should match 12a)
#   --great_background          all (grey included; default) | assigned (non-grey only)
#   --min_gene_set_size / --max_gene_set_size   gene-set size filter [5 / Inf for GREAT; cP uses 10-500]
#   --great_engine              fast (default; rGREAT statistics, coverage computed once) | rgreat
#   --min_region_hits / --min_gene_hits         GREAT: terms FDR-tested if >= 5 region hits (rGREAT
#                               default); flagged significant if FDR < --fdr_cut and >= 2 genes
#   --fdr_cut                   FDR for significance flags and plots [0.05]
#   --enrichr_dbs, --use_enrichr_background [TRUE], --min_bg_genes, --enrichr_retries, --enrichr_sleep
#   --plot_p_col, --plot_top_n, --plot_sig_cutoff, --plot_width, --plot_height
#   --do_kegg / --do_go / --do_reactome   (kept for old SLURM scripts: FALSE drops those collections)
#   --organism                  (kept for old SLURM scripts; ignored -- KEGG comes from the cache)
#
# Outputs:
#   <project_root>/comethyl_output/12b_enrichment/<cpg>/<region>/<variant>/[<list_stem>/]
#     selected_modules.txt, selected_modules_summary.tsv, enrichment_overview.tsv
#     GREAT/            <module>_GREAT_<collection>.tsv/.xlsx, Summary_GREAT_<collection>_*.pdf/.png
#     clusterProfiler/  <module>_<collection>.tsv/.xlsx,       Summary_<collection>_*.pdf/.png
#     Enrichr/          <module>_enrichr.xlsx, <module>_enrichr_<DB>.pdf, Summary_Enrichr_<DB>_*.pdf/.png,
#                       enrichr_status.tsv
#     run_log.txt, run_parameters.txt, session_info.txt
# ============================================================

message("Starting Script 12b")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(scales)
  library(GenomicRanges)
})

# ============================================================
# 1) CLI
# ============================================================
option_list <- list(
  make_option("--project_root", type = "character", default = NULL,
              help = "Project root. Outputs go under <project_root>/comethyl_output/12b_enrichment/. Required unless --build_gene_set_cache_only TRUE."),
  make_option("--annotation_dir_v1", type = "character", default = NULL,
              help = "Script 12a output folder for variant v1 (contains Annotated_Regions.tsv, Module_Gene_List.tsv, Background_Genes.tsv)."),
  make_option("--annotation_dir_v2", type = "character", default = NULL, help = "12a output folder for variant v2."),
  make_option("--annotation_dir_v3", type = "character", default = NULL, help = "12a output folder for variant v3."),
  make_option("--modules", type = "character", default = NULL,
              help = "Comma-separated modules to enrich, e.g. 'yellow,brown,red,lightsteelblue1'. Default (with no list/trait file): all non-grey modules."),
  make_option("--module_list_file", type = "character", default = NULL,
              help = "Text file of module names, or a directory of .txt files (each processed into a subfolder named after the file)."),
  make_option("--trait_list_file", type = "character", default = NULL,
              help = "Text file (or directory) of traits; modules with ME-trait p <= --p_thresh in the stats file(s) are selected."),
  make_option("--stats_file_v1", type = "character", default = NULL, help = "ME-trait stats (Script 09a) for v1."),
  make_option("--stats_file_v2", type = "character", default = NULL, help = "ME-trait stats for v2."),
  make_option("--stats_file_v3", type = "character", default = NULL, help = "ME-trait stats for v3."),
  make_option("--p_thresh", type = "double", default = 0.05, help = "Module selection p threshold for trait mode [default %default]."),
  make_option("--collections", type = "character", default = "all",
              help = "Gene-set collections for GREAT and clusterProfiler: comma list of GO_BP,GO_MF,GO_CC,KEGG,Reactome,Hallmark, or 'all' [default %default]."),
  make_option("--gene_set_cache", type = "character", default = NULL,
              help = "Gene-set cache RDS (Entrez gene sets for all collections). Build once with --build_gene_set_cache_only TRUE."),
  make_option("--build_gene_set_cache_only", type = "character", default = "FALSE",
              help = "TRUE to build --gene_set_cache and exit. Needs internet for KEGG/MSigDB (run on a login node) [default %default]."),
  make_option("--gene_set_source", type = "character", default = "native",
              help = "Cache source: 'native' (GO org.Hs.eg.db, KEGG REST, reactome.db, Hallmark msigdbr) or 'msigdb' (all from one MSigDB release) [default %default]."),
  make_option("--do_great", type = "character", default = "TRUE", help = "Run GREAT region-based enrichment [default %default]."),
  make_option("--do_clusterprofiler", type = "character", default = "TRUE",
              help = "Run clusterProfiler ORA with universe = Background_Genes [default %default]."),
  make_option("--do_enrichr", type = "character", default = "TRUE", help = "Run Enrichr (needs internet) [default %default]."),
  make_option("--tss_source", type = "character", default = "GREAT:hg38",
              help = "GREAT gene/TSS set: 'GREAT:hg38' (GREAT v4 genes; comethyl default) or 'TxDb:hg38'. Should match 12a [default %default]."),
  make_option("--great_mode", type = "character", default = "basalPlusExt",
              help = "GREAT association rule: basalPlusExt, twoClosest or oneClosest [default %default]."),
  make_option("--basal_upstream", type = "double", default = 5000, help = "basalPlusExt upstream bp [default %default]."),
  make_option("--basal_downstream", type = "double", default = 1000, help = "basalPlusExt downstream bp [default %default]."),
  make_option("--extension", type = "double", default = 1e6, help = "Maximum domain extension bp [default %default]."),
  make_option("--great_background", type = "character", default = "all",
              help = "GREAT background regions: 'all' = every region incl. grey; 'assigned' = non-grey only [default %default]."),
  make_option("--min_gene_set_size", type = "double", default = 5, help = "Minimum gene-set size for GREAT [default %default]."),
  make_option("--max_gene_set_size", type = "double", default = Inf, help = "Maximum gene-set size for GREAT [default %default]."),
  make_option("--great_engine", type = "character", default = "fast",
              help = "'fast' = rGREAT's statistics with domain coverage computed once for all modules (seconds); 'rgreat' = call rGREAT::great() per module (minutes per module; for verification) [default %default]."),
  make_option("--min_region_hits", type = "double", default = 5,
              help = "GREAT: terms need >= this many module region hits to be FDR-tested and flagged (rGREAT getEnrichmentTable() default) [default %default]."),
  make_option("--min_gene_hits", type = "double", default = 2, help = "GREAT: minimum term genes hit to flag a term [default %default]."),
  make_option("--fdr_cut", type = "double", default = 0.05, help = "FDR threshold for significance flags [default %default]."),
  make_option("--enrichr_dbs", type = "character",
              default = "GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,Reactome_Pathways_2024,MSigDB_Hallmark_2020,GWAS_Catalog_2025",
              help = "Comma-separated Enrichr libraries; validated against Enrichr's current list [default %default]."),
  make_option("--use_enrichr_background", type = "character", default = "TRUE",
              help = "Use 12a's Background_Genes as the Enrichr background [default %default]."),
  make_option("--min_bg_genes", type = "integer", default = 200L, help = "Minimum background size to use it [default %default]."),
  make_option("--enrichr_retries", type = "integer", default = 3L, help = "Enrichr retries per module [default %default]."),
  make_option("--enrichr_sleep", type = "double", default = 1, help = "Seconds between Enrichr calls [default %default]."),
  make_option("--plot_p_col", type = "character", default = "p.adjust", help = "'p.adjust' or 'pvalue' for summary plots [default %default]."),
  make_option("--plot_top_n", type = "integer", default = 10L, help = "Top N terms per module in summary plots [default %default]."),
  make_option("--plot_sig_cutoff", type = "double", default = 0.05, help = "Significance cutoff for summary plots [default %default]."),
  make_option("--plot_width", type = "double", default = 14, help = "Summary plot width, inches [default %default]."),
  make_option("--plot_height", type = "double", default = 10, help = "Summary plot height, inches [default %default]."),
  make_option("--helper_file", type = "character", default = NULL, help = "Path to helper.R [default: next to this script]."),
  make_option("--do_kegg", type = "character", default = NULL, help = "Deprecated: FALSE drops KEGG from --collections."),
  make_option("--do_go", type = "character", default = NULL, help = "Deprecated: FALSE drops GO_BP/MF/CC from --collections."),
  make_option("--do_reactome", type = "character", default = NULL, help = "Deprecated: FALSE drops Reactome from --collections."),
  make_option("--organism", type = "character", default = NULL, help = "Deprecated; ignored (KEGG comes from the gene-set cache).")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "GREAT, clusterProfiler and Enrichr enrichment of comethyl modules."))

to_bool <- function(x, nm) {
  v <- toupper(trimws(as.character(x)))
  if (v %in% c("TRUE", "T", "1", "YES", "Y")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO", "N")) return(FALSE)
  stop("--", nm, " must be TRUE or FALSE (got '", x, "')", call. = FALSE)
}
for (nm in c("build_gene_set_cache_only", "do_great", "do_clusterprofiler", "do_enrichr", "use_enrichr_background")) {
  opt[[nm]] <- to_bool(opt[[nm]], nm)
}

# ============================================================
# 2) Helper
# ============================================================
script_dir <- tryCatch({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
helper_file <- if (!is.null(opt$helper_file)) opt$helper_file else file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found at: ", helper_file, call. = FALSE)
source(helper_file)
if (exists("setup_publication_font", mode = "function")) setup_publication_font()

# ============================================================
# 3) Collections / validation
# ============================================================
collections <- parse_collections(opt$collections)
if (!is.null(opt$do_kegg) && !to_bool(opt$do_kegg, "do_kegg")) collections <- setdiff(collections, "KEGG")
if (!is.null(opt$do_go) && !to_bool(opt$do_go, "do_go")) collections <- setdiff(collections, c("GO_BP", "GO_MF", "GO_CC"))
if (!is.null(opt$do_reactome) && !to_bool(opt$do_reactome, "do_reactome")) collections <- setdiff(collections, "Reactome")
if (!is.null(opt$organism)) message("Note: --organism is ignored (KEGG gene sets come from the gene-set cache).")
opt$gene_set_source <- match.arg(opt$gene_set_source, c("native", "msigdb"))
opt$great_mode <- match.arg(opt$great_mode, c("basalPlusExt", "twoClosest", "oneClosest"))
opt$great_background <- match.arg(opt$great_background, c("all", "assigned"))
opt$great_engine <- match.arg(opt$great_engine, c("fast", "rgreat"))
if (!opt$plot_p_col %in% c("p.adjust", "pvalue")) stop("--plot_p_col must be 'p.adjust' or 'pvalue'", call. = FALSE)

if (opt$build_gene_set_cache_only) {
  stop_if_missing(opt$gene_set_cache, "--gene_set_cache")
  build_gene_set_cache(collections, source = opt$gene_set_source, cache_file = opt$gene_set_cache)
  message("Script 12b: gene-set cache built: ", opt$gene_set_cache)
  quit(save = "no", status = 0)
}

stop_if_missing(opt$project_root, "--project_root")
stop_if_missing(opt$annotation_dir_v1, "--annotation_dir_v1")
if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root, call. = FALSE)
for (nm in c("annotation_dir_v1", "annotation_dir_v2", "annotation_dir_v3")) {
  if (!is.null(opt[[nm]]) && !dir.exists(opt[[nm]])) stop("--", nm, " not found: ", opt[[nm]], call. = FALSE)
}
need_sets <- opt$do_great || opt$do_clusterprofiler
gene_sets <- NULL
if (need_sets) {
  stop_if_missing(opt$gene_set_cache, "--gene_set_cache (needed for GREAT / clusterProfiler)")
  gene_sets <- load_gene_set_cache(opt$gene_set_cache, collections)
}
if (length(collections) == 0 && need_sets) stop("No collections left to test.", call. = FALSE)

# ============================================================
# 4) Module-selection helpers (as in the previous 12b)
# ============================================================
resolve_list_input <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (file.exists(path) && !dir.exists(path)) return(stats::setNames(normalizePath(path), ""))
  if (dir.exists(path)) {
    txts <- sort(list.files(path, pattern = "\\.txt$", full.names = TRUE))
    if (length(txts) == 0) stop("Directory contains no .txt files: ", path, call. = FALSE)
    return(stats::setNames(normalizePath(txts), tools::file_path_sans_ext(basename(txts))))
  }
  stop("Path does not exist: ", path, call. = FALSE)
}
read_clean_txt_list <- function(path) {
  x <- trimws(readLines(path, warn = FALSE)); x <- x[nzchar(x) & !grepl("^#", x)]; unique(x)
}
read_stats_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext, xlsx = openxlsx::read.xlsx(path), xls = openxlsx::read.xlsx(path),
               csv = utils::read.csv(path, check.names = FALSE),
               utils::read.delim(path, check.names = FALSE))
  cn <- tolower(colnames(df))
  for (col in c("module", "trait")) {
    if (!col %in% cn) stop("stats file must contain a ", col, " column: ", path, call. = FALSE)
    names(df)[match(col, cn)] <- col
  }
  pc <- c("p", "pvalue", "p.value", "p_value")[c("p", "pvalue", "p.value", "p_value") %in% tolower(colnames(df))][1]
  if (is.na(pc)) stop("stats file needs a p-value column (p/pvalue/p.value/p_value): ", path, call. = FALSE)
  names(df)[match(pc, tolower(colnames(df)))] <- "p"
  df
}

load_annotation_files <- function(annotation_dir) {
  f <- function(x) file.path(annotation_dir, x)
  for (x in c("Annotated_Regions.tsv", "Module_Gene_List.tsv", "Module_Gene_Summary.tsv", "Background_Genes.tsv")) {
    validate_file_exists(f(x), x)
  }
  regions <- utils::read.delim(f("Annotated_Regions.tsv"), stringsAsFactors = FALSE, check.names = FALSE,
                               colClasses = c(chr = "character", module = "character"))
  miss <- setdiff(c("RegionID", "chr", "start", "end", "module"), colnames(regions))
  if (length(miss) > 0) stop("Annotated_Regions.tsv missing: ", paste(miss, collapse = ", "), call. = FALSE)
  params <- NULL
  if (file.exists(f("run_parameters.txt"))) {
    p <- utils::read.delim(f("run_parameters.txt"), header = FALSE, stringsAsFactors = FALSE, quote = "")
    params <- stats::setNames(as.character(p[[2]]), p[[1]])
  }
  list(regions = regions[, c("RegionID", "chr", "start", "end", "module")],
       gene_list = utils::read.delim(f("Module_Gene_List.tsv"), stringsAsFactors = FALSE),
       module_summary = utils::read.delim(f("Module_Gene_Summary.tsv"), stringsAsFactors = FALSE),
       background_genes = utils::read.delim(f("Background_Genes.tsv"), stringsAsFactors = FALSE)$gene_symbol,
       params = params)
}

# ============================================================
# 5) Plot helpers (previous 12b dot-plot style; PDF via cairo + PNG)
# ============================================================
make_enrichment_dotplot <- function(df, title, source_tag = NULL, p_col_label = "adj. p",
                                    size_label = "Gene Count", sig_cutoff = 0.05, top_n = 10,
                                    max_chars = 55, plot_width = 14, plot_height = 10,
                                    base_size = 13, file_stem = NULL) {
  if (is.null(df) || nrow(df) == 0) { message("[plot] No data for: ", title); return(invisible(NULL)) }
  df$Description <- ifelse(nchar(df$Description) > max_chars,
                           paste0(substr(df$Description, 1, max_chars - 3), "..."), df$Description)
  df <- df %>% dplyr::filter(!is.na(p_col_val), p_col_val <= sig_cutoff) %>%
    dplyr::group_by(module) %>% dplyr::arrange(p_col_val, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>% dplyr::ungroup()
  if (nrow(df) == 0) { message("[plot] No terms pass ", sig_cutoff, " for: ", title); return(invisible(NULL)) }
  df <- df %>% dplyr::arrange(module, p_col_val) %>%
    dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))
  p <- ggplot(df, aes(x = module, y = Description, size = gene_count, fill = p_col_val)) +
    geom_point(shape = 23, color = "grey30") +
    scale_size(range = c(3, 10), name = size_label) +
    scale_fill_gradient(low = "#B22222", high = "#FFC0CB", name = p_col_label,
                        limits = c(0, sig_cutoff), oob = scales::squish) +
    theme_bw(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = base_size - 1),
          axis.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = title, subtitle = paste0("Top ", top_n, " terms per module | ", p_col_label, " <= ", sig_cutoff,
                                          if (!is.null(source_tag)) paste0(" | Source: ", source_tag) else ""))
  if (!is.null(file_stem)) save_annotation_plot(p, file_stem, plot_width, plot_height)
  invisible(p)
}

write_tsv_xlsx <- function(df, stem) {
  utils::write.table(df, paste0(stem, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  openxlsx::write.xlsx(df, paste0(stem, ".xlsx"), rowNames = FALSE)
}

# ============================================================
# 6) GREAT (region-based, offline)
# ============================================================
extended_tss <- NULL
if (opt$do_great) {
  extended_tss <- great_build_extended_tss(opt$tss_source, genome = "hg38", mode = opt$great_mode,
                                           basal_upstream = opt$basal_upstream,
                                           basal_downstream = opt$basal_downstream,
                                           extension = opt$extension, verbose = TRUE)
}

# All selected collections combined into one gene-set list (one great() call
# per module), names "<collection>::<term_id>", size-filtered.
combined_gene_sets <- function(gs, collections, min_size, max_size, universe = NULL) {
  out <- list()
  for (coll in collections) {
    s <- gs$collections[[coll]]
    if (!is.null(universe)) s <- lapply(s, intersect, universe)
    n <- lengths(s)
    s <- s[n >= min_size & n <= max_size]
    if (length(s) > 0) names(s) <- paste0(coll, "::", names(s))
    out <- c(out, s)
  }
  out
}

# rGREAT::great() per module (verification engine). Returns rGREAT's full
# table (no min-hit filter) plus the genes hit by the module.
run_great_rgreat <- function(module_regions, background_gr, great_sets, ext, min_size) {
  gr <- GenomicRanges::GRanges(module_regions$chr, IRanges::IRanges(module_regions$start, module_regions$end))
  res <- rGREAT::great(gr, gene_sets = great_sets, extended_tss = ext, background = background_gr,
                       exclude = NULL, min_gene_set_size = min_size, verbose = FALSE)
  tb <- rGREAT::getEnrichmentTable(res, min_region_hits = 0)
  mid <- GenomicRanges::GRanges(module_regions$chr,
                                IRanges::IRanges(floor((module_regions$start + module_regions$end) / 2), width = 1))
  attr(tb, "hit_genes") <- unique(ext$gene_id[S4Vectors::subjectHits(GenomicRanges::findOverlaps(mid, ext))])
  tb
}

# Add collection / term name / genes, FDR within module x collection over
# terms with >= min_region_hits (rGREAT getEnrichmentTable() convention).
finish_great_table <- function(tb, great_sets, ext, gs, coll) {
  hit_genes <- attr(tb, "hit_genes"); sym <- attr(ext, "symbol_map")
  tb$collection <- sub("::.*$", "", tb$id)
  d <- tb[tb$collection == coll, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d$term_id <- sub("^[^:]+::", "", d$id)
  nm <- gs$term_names[[coll]][d$term_id]
  d$Description <- ifelse(is.na(nm), d$term_id, unname(nm))
  d$genes <- vapply(d$id, function(i) {
    g <- intersect(great_sets[[i]], hit_genes); s2 <- sym[g]; s2[is.na(s2)] <- g[is.na(s2)]
    paste(sort(unname(s2)), collapse = ",")
  }, "")
  tested <- d$observed_region_hits >= opt$min_region_hits
  d$fdr <- NA_real_; d$fdr_hyper <- NA_real_
  d$fdr[tested] <- stats::p.adjust(d$p_value[tested], "BH")
  d$fdr_hyper[tested] <- stats::p.adjust(d$p_value_hyper[tested], "BH")
  d$significant <- tested & !is.na(d$fdr) & d$fdr < opt$fdr_cut & d$observed_gene_hits >= opt$min_gene_hits
  d <- d[order(is.na(d$fdr), d$fdr, d$p_value),
         c("collection", "term_id", "Description", "genome_fraction", "observed_region_hits", "fold_enrichment",
           "p_value", "fdr", "observed_gene_hits", "gene_set_size", "fold_enrichment_hyper", "p_value_hyper",
           "fdr_hyper", "mean_tss_dist", "significant", "genes")]
  rownames(d) <- NULL
  d
}

# ============================================================
# 7) clusterProfiler ORA with universe
# ============================================================
run_clusterprofiler_module <- function(gene_symbols, universe_entrez, gs, collections, out_prefix, log_fn) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) { log_fn("clusterProfiler not installed; skipping."); return(invisible(NULL)) }
  fg <- symbols_to_entrez(gene_symbols, gs)
  fg_e <- intersect(fg$entrez, universe_entrez)
  log_fn("  clusterProfiler: ", fg$n_mapped, "/", fg$n_in, " genes mapped to Entrez; ", length(fg_e), " in universe")
  if (length(fg_e) == 0) return(invisible(NULL))
  eg2sym <- stats::setNames(names(gs$symbol_map), gs$symbol_map)
  for (coll in collections) {
    s <- gs$collections[[coll]]
    t2g <- data.frame(term = rep(names(s), lengths(s)), gene = unlist(s, use.names = FALSE), stringsAsFactors = FALSE)
    t2n <- data.frame(term = names(gs$term_names[[coll]]), name = unname(gs$term_names[[coll]]), stringsAsFactors = FALSE)
    res <- tryCatch(suppressMessages(clusterProfiler::enricher(fg_e, TERM2GENE = t2g, TERM2NAME = t2n,
                                                               universe = universe_entrez,
                                                               pvalueCutoff = 1, qvalueCutoff = 1)),
                    error = function(e) { log_fn("  ERROR clusterProfiler ", coll, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) next
    df <- as.data.frame(res)
    if (nrow(df) == 0) next
    df$geneSymbols <- vapply(strsplit(df$geneID, "/"), function(g) {
      s2 <- eg2sym[g]; s2[is.na(s2)] <- g[is.na(s2)]; paste(unname(s2), collapse = "/")
    }, "")
    df$collection <- coll
    write_tsv_xlsx(df, paste0(out_prefix, "_", coll))
  }
  invisible(NULL)
}

# ============================================================
# 8) Enrichr (online)
# ============================================================
init_enrichr <- function() {
  options(enrichR.sites.base.address = "https://maayanlab.cloud/",
          enrichR.base.address = "https://maayanlab.cloud/Enrichr/",
          speedrichr.base.address = "https://maayanlab.cloud/speedrichr/api/",
          enrichR.live = TRUE, enrichR.quiet = TRUE, modEnrichR.use = TRUE,
          enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr", "OxEnrichr"))
}
# Thin wrappers so the API calls live in one place (defined only if not
# already present, so a test harness sourced via --helper_file can mock them).
if (!exists("enrichr_list_dbs", mode = "function")) {
  enrichr_list_dbs <- function() { init_enrichr(); enrichR::listEnrichrDbs() }
}
if (!exists("enrichr_query", mode = "function")) {
  enrichr_query <- function(genes, dbs, background = NULL) {
    init_enrichr()
    if (is.null(background)) enrichR::enrichr(genes = genes, databases = dbs)
    else enrichR::enrichr(genes = genes, databases = dbs, background = background, include_overlap = TRUE)
  }
}

enrichr_preflight <- function(dbs, want_background) {
  out <- list(ok = FALSE, dbs = character(0), background_supported = FALSE, reason = "")
  if (!requireNamespace("enrichR", quietly = TRUE)) { out$reason <- "enrichR package not installed"; return(out) }
  ver <- utils::packageVersion("enrichR")
  out$background_supported <- ver >= "3.3" && "background" %in% names(formals(enrichR::enrichr))
  if (want_background && !out$background_supported) {
    out$reason <- paste0("enrichR ", ver, " has no background support (need >= 3.3)")
    return(out)
  }
  avail <- tryCatch(enrichr_list_dbs(), error = function(e) e)
  if (inherits(avail, "error") || is.null(avail) || !"libraryName" %in% colnames(avail)) {
    out$reason <- paste0("Enrichr unreachable: ",
                         if (inherits(avail, "error")) conditionMessage(avail) else "no library list returned")
    return(out)
  }
  out$dbs <- intersect(dbs, avail$libraryName)
  out$dropped <- setdiff(dbs, avail$libraryName)
  out$ok <- length(out$dbs) > 0
  if (!out$ok) out$reason <- "none of the requested Enrichr libraries exist"
  out
}

run_enrichr_module <- function(gene_symbols, background, dbs, out_prefix, mod_name, retries, sleep_time, log_fn) {
  genes <- unique(trimws(gene_symbols)); genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) return(list(status = "skipped", detail = "no genes"))
  res <- NULL; last_err <- NULL
  for (attempt in seq_len(retries)) {
    res <- tryCatch(enrichr_query(genes, dbs, background), error = function(e) { last_err <<- conditionMessage(e); NULL })
    if (!is.null(res) && length(res) > 0) break
    if (attempt < retries) Sys.sleep(sleep_time)
  }
  if (is.null(res) || length(res) == 0) {
    log_fn("  ERROR Enrichr failed for ", mod_name, ": ", if (is.null(last_err)) "empty result" else last_err)
    return(list(status = "failed", detail = if (is.null(last_err)) "empty result" else last_err))
  }
  wb <- openxlsx::createWorkbook(); n_rows <- 0
  for (nm in names(res)) {
    df <- res[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    if (!"Overlap" %in% colnames(df)) {  # Speedrichr without overlap: derive count from Genes
      df$Overlap <- paste0(lengths(strsplit(as.character(df$Genes), ";")), "/NA")
    }
    sheet <- substr(gsub("[^A-Za-z0-9]", "_", nm), 1, 31)
    openxlsx::addWorksheet(wb, sheet); openxlsx::writeData(wb, sheet, df)
    n_rows <- n_rows + nrow(df)
    pdf_stem <- paste0(out_prefix, "_enrichr_", substr(gsub("[^A-Za-z0-9]", "_", nm), 1, 60))
    tryCatch({
      p <- enrichR::plotEnrich(df, showTerms = 25, numChar = 75, y = "Count", orderBy = "P.value") +
        ggplot2::ggtitle(paste(nm, "for", mod_name, "module"))
      save_annotation_plot(p, pdf_stem, width = 15, height = 7, verbose = FALSE, png = FALSE)
    }, error = function(e) log_fn("  WARNING Enrichr plot failed for ", nm, ": ", conditionMessage(e)))
  }
  if (n_rows > 0) openxlsx::saveWorkbook(wb, paste0(out_prefix, "_enrichr.xlsx"), overwrite = TRUE)
  Sys.sleep(sleep_time)
  list(status = "ok", detail = paste0(n_rows, " rows across ", length(res), " libraries"))
}

# ============================================================
# 9) Per-variant runner
# ============================================================
run_enrichment_for_variant <- function(annotation_dir, variant_label, stats_file = NULL,
                                       module_list_file = NULL, trait_list = NULL, list_stem = "",
                                       enrichr_state = NULL) {
  files <- load_annotation_files(annotation_dir)
  variant_name <- basename(annotation_dir); region_label <- basename(dirname(annotation_dir))
  cpg_label <- basename(dirname(dirname(annotation_dir)))
  out_dir <- file.path(opt$project_root, "comethyl_output", "12b_enrichment", cpg_label, region_label, variant_name)
  if (nzchar(list_stem)) out_dir <- file.path(out_dir, list_stem)
  for (d in c("GREAT", "clusterProfiler", "Enrichr")) safe_dir_create(file.path(out_dir, d))
  log_file <- file.path(out_dir, "run_log.txt")
  L <- function(...) append_log(log_file, ...)
  L("Starting enrichment for ", variant_label, if (nzchar(list_stem)) paste0(" / list ", list_stem) else "")
  L("annotation_dir: ", annotation_dir)
  L("collections: ", paste(collections, collapse = ", "), "; GREAT: ", opt$do_great,
    "; clusterProfiler: ", opt$do_clusterprofiler, "; Enrichr: ", opt$do_enrichr)

  # 12a/12b GREAT settings consistency
  pa <- files$params
  if (!is.null(pa) && "great_tss_source" %in% names(pa)) {
    cmp <- c(tss_source = opt$tss_source, mode = opt$great_mode)
    got <- c(tss_source = pa[["great_tss_source"]], mode = pa[["great_mode"]])
    if (!identical(tolower(unname(cmp)), tolower(unname(got)))) {
      L("WARNING: GREAT settings differ from 12a (12a: ", paste(got, collapse = "/"), "; 12b: ",
        paste(cmp, collapse = "/"), "). Term genes will not match 12a's gene columns.")
    }
  } else if (!is.null(pa) && "gene_assignment" %in% names(pa) && pa[["gene_assignment"]] != "great") {
    L("WARNING: 12a used gene_assignment = ", pa[["gene_assignment"]], "; Enrichr/clusterProfiler gene lists are not GREAT-based.")
  }

  # module selection
  all_mods <- sort(setdiff(unique(files$regions$module), "grey"))
  if (!is.null(module_list_file)) {
    sel <- read_clean_txt_list(module_list_file); mode_sel <- "module_list_file"
  } else if (!is.null(trait_list)) {
    st <- read_stats_file(stats_file)
    sel <- unique(as.character(st$module[!is.na(st$p) & st$trait %in% trait_list & as.numeric(st$p) <= opt$p_thresh]))
    mode_sel <- "stats_file_plus_trait_list"
  } else if (!is.null(opt$modules)) {
    sel <- trimws(strsplit(opt$modules, ",")[[1]]); mode_sel <- "--modules"
  } else {
    sel <- all_mods; mode_sel <- "all_modules"
  }
  missing_mods <- setdiff(sel, all_mods)
  if (length(missing_mods) > 0) L("WARNING: modules not found (or grey) and skipped: ", paste(missing_mods, collapse = ", "))
  sel <- intersect(sel, all_mods)
  write_lines_safe(sel, file.path(out_dir, "selected_modules.txt"))
  utils::write.table(files$module_summary[files$module_summary$module %in% sel, ],
                     file.path(out_dir, "selected_modules_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  L("selection: ", mode_sel, "; ", length(sel), " modules: ", paste(sel, collapse = ", "))
  if (length(sel) == 0) { L("No modules selected; nothing to do."); return(invisible(NULL)) }

  overview <- list()

  # ---- GREAT ----
  if (opt$do_great) {
    bg_df <- if (opt$great_background == "all") files$regions else files$regions[files$regions$module != "grey", ]
    bg_gr <- GenomicRanges::GRanges(bg_df$chr, IRanges::IRanges(bg_df$start, bg_df$end))
    great_sets <- combined_gene_sets(gene_sets, collections, opt$min_gene_set_size, opt$max_gene_set_size,
                                     universe = unique(extended_tss$gene_id))
    L("GREAT: background = ", opt$great_background, " (", nrow(bg_df), " regions); ", length(great_sets),
      " gene sets after size filter (", opt$min_gene_set_size, "-", opt$max_gene_set_size, ", genes in TSS set)")
    t0 <- Sys.time()
    tables <- list()
    if (opt$great_engine == "fast") {
      prep <- great_fast_prepare(extended_tss, bg_gr, great_sets, opt$min_gene_set_size, opt$max_gene_set_size)
      tables <- great_fast_test(prep, split(files$regions[files$regions$module %in% sel, ],
                                            files$regions$module[files$regions$module %in% sel]))
    } else {
      for (m in sel) {
        tables[[m]] <- tryCatch(run_great_rgreat(files$regions[files$regions$module == m, ], bg_gr, great_sets,
                                                 extended_tss, opt$min_gene_set_size),
                                error = function(e) { L("  ERROR GREAT ", m, ": ", conditionMessage(e)); NULL })
      }
    }
    L("GREAT (", opt$great_engine, " engine): ", length(sel), " modules in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
    for (m in names(tables)) {
      if (is.null(tables[[m]])) next
      for (coll in collections) {
        d <- finish_great_table(tables[[m]], great_sets, extended_tss, gene_sets, coll)
        if (is.null(d)) next
        write_tsv_xlsx(d, file.path(out_dir, "GREAT", paste0(m, "_GREAT_", coll)))
        overview[[length(overview) + 1]] <- data.frame(module = m, method = "GREAT", collection = coll,
                                                       n_tested = sum(!is.na(d$fdr)), n_significant = sum(d$significant))
      }
    }
    for (coll in collections) {
      dd <- do.call(rbind, lapply(sel, function(m) {
        f <- file.path(out_dir, "GREAT", paste0(m, "_GREAT_", coll, ".tsv"))
        if (!file.exists(f)) return(NULL)
        d <- utils::read.delim(f, stringsAsFactors = FALSE, quote = "")
        d <- d[!is.na(d$fdr) & d$observed_gene_hits >= opt$min_gene_hits, ]
        if (nrow(d) == 0) return(NULL)
        data.frame(module = m, Description = d$Description,
                   p_col_val = if (opt$plot_p_col == "p.adjust") d$fdr else d$p_value,
                   gene_count = d$observed_region_hits)
      }))
      make_enrichment_dotplot(dd, paste0("Top GREAT ", coll, " terms per module"), source_tag = "rGREAT (local)",
                              p_col_label = if (opt$plot_p_col == "p.adjust") "FDR" else "raw p",
                              size_label = "Region hits", sig_cutoff = opt$plot_sig_cutoff, top_n = opt$plot_top_n,
                              plot_width = opt$plot_width, plot_height = opt$plot_height,
                              file_stem = file.path(out_dir, "GREAT", sprintf("Summary_GREAT_%s_top%d_%s%s", coll,
                                opt$plot_top_n, gsub("\\.", "", opt$plot_p_col), opt$plot_sig_cutoff)))
    }
  }

  # ---- background genes (shared by clusterProfiler and Enrichr) ----
  bg_syms <- unique(files$background_genes[!is.na(files$background_genes) & nzchar(files$background_genes)])
  L("Background genes (12a, genes assigned to any region): ", length(bg_syms))

  # ---- clusterProfiler ----
  if (opt$do_clusterprofiler) {
    uni <- symbols_to_entrez(bg_syms, gene_sets)
    L("clusterProfiler universe: ", uni$n_mapped, "/", uni$n_in, " background genes mapped to Entrez")
    for (m in sel) {
      genes <- files$gene_list$gene_symbol[files$gene_list$module == m]
      L("Module ", m, ": ", length(unique(genes)), " genes")
      run_clusterprofiler_module(genes, uni$entrez, gene_sets, collections,
                                 file.path(out_dir, "clusterProfiler", m), L)
      for (coll in collections) {
        f <- file.path(out_dir, "clusterProfiler", paste0(m, "_", coll, ".tsv"))
        n_sig <- if (file.exists(f)) sum(utils::read.delim(f, quote = "")$p.adjust < opt$fdr_cut, na.rm = TRUE) else NA
        overview[[length(overview) + 1]] <- data.frame(module = m, method = "clusterProfiler", collection = coll,
                                                       n_tested = if (file.exists(f)) nrow(utils::read.delim(f, quote = "")) else 0,
                                                       n_significant = n_sig)
      }
    }
    for (coll in collections) {
      dd <- do.call(rbind, lapply(sel, function(m) {
        f <- file.path(out_dir, "clusterProfiler", paste0(m, "_", coll, ".tsv"))
        if (!file.exists(f)) return(NULL)
        d <- utils::read.delim(f, stringsAsFactors = FALSE, quote = "")
        data.frame(module = m, Description = d$Description,
                   p_col_val = if (opt$plot_p_col == "p.adjust") d$p.adjust else d$pvalue, gene_count = d$Count)
      }))
      make_enrichment_dotplot(dd, paste0("Top ", coll, " terms per module"), source_tag = "clusterProfiler (universe = background)",
                              p_col_label = if (opt$plot_p_col == "p.adjust") "adj. p" else "raw p",
                              sig_cutoff = opt$plot_sig_cutoff, top_n = opt$plot_top_n,
                              plot_width = opt$plot_width, plot_height = opt$plot_height,
                              file_stem = file.path(out_dir, "clusterProfiler", sprintf("Summary_%s_top%d_%s%s", coll,
                                opt$plot_top_n, gsub("\\.", "", opt$plot_p_col), opt$plot_sig_cutoff)))
    }
  }

  # ---- Enrichr ----
  if (opt$do_enrichr) {
    status <- list()
    if (!isTRUE(enrichr_state$ok)) {
      L("ERROR: Enrichr skipped -- ", enrichr_state$reason)
      status[[1]] <- data.frame(module = sel, status = "not_run", detail = enrichr_state$reason)
    } else {
      bg_use <- NULL
      if (opt$use_enrichr_background) {
        if (length(bg_syms) >= opt$min_bg_genes) bg_use <- bg_syms
        else L("WARNING: background has ", length(bg_syms), " genes (< ", opt$min_bg_genes, "); Enrichr runs WITHOUT background")
      }
      L("Enrichr: libraries = ", paste(enrichr_state$dbs, collapse = ", "),
        "; background = ", if (is.null(bg_use)) "none" else paste0(length(bg_use), " genes"))
      for (m in sel) {
        genes <- files$gene_list$gene_symbol[files$gene_list$module == m]
        r <- run_enrichr_module(genes, bg_use, enrichr_state$dbs, file.path(out_dir, "Enrichr", m), m,
                                opt$enrichr_retries, opt$enrichr_sleep, L)
        status[[length(status) + 1]] <- data.frame(module = m, status = r$status, detail = r$detail,
                                                   background_used = !is.null(bg_use))
      }
      for (db in enrichr_state$dbs) {
        sheet <- substr(gsub("[^A-Za-z0-9]", "_", db), 1, 31)
        dd <- do.call(rbind, lapply(sel, function(m) {
          f <- file.path(out_dir, "Enrichr", paste0(m, "_enrichr.xlsx"))
          if (!file.exists(f) || !sheet %in% openxlsx::getSheetNames(f)) return(NULL)
          d <- openxlsx::read.xlsx(f, sheet = sheet)
          pc <- if (opt$plot_p_col == "p.adjust") "Adjusted.P.value" else "P.value"
          data.frame(module = m, Description = gsub(" \\(GO:\\d+\\)$", "", d$Term),
                     p_col_val = as.numeric(d[[pc]]),
                     gene_count = suppressWarnings(as.integer(sub("/.*$", "", as.character(d$Overlap)))))
        }))
        make_enrichment_dotplot(dd, paste0("Top Enrichr terms: ", db), source_tag = "Enrichr / maayanlab.cloud",
                                p_col_label = if (opt$plot_p_col == "p.adjust") "adj. p" else "raw p",
                                sig_cutoff = opt$plot_sig_cutoff, top_n = opt$plot_top_n,
                                plot_width = opt$plot_width, plot_height = opt$plot_height,
                                file_stem = file.path(out_dir, "Enrichr", sprintf("Summary_Enrichr_%s_top%d_%s%s",
                                  gsub("[^A-Za-z0-9]+", "_", db), opt$plot_top_n, gsub("\\.", "", opt$plot_p_col), opt$plot_sig_cutoff)))
      }
    }
    st <- do.call(rbind, status)
    utils::write.table(st, file.path(out_dir, "Enrichr", "enrichr_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    n_bad <- sum(st$status != "ok")
    if (n_bad > 0) L("ERROR: Enrichr did not complete for ", n_bad, " of ", nrow(st), " modules -- see Enrichr/enrichr_status.tsv")
  }

  if (length(overview) > 0) {
    utils::write.table(do.call(rbind, overview), file.path(out_dir, "enrichment_overview.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }

  params <- c(paste0("timestamp\t", timestamp_now()), paste0("variant_label\t", variant_label),
              paste0("list_stem\t", list_stem), paste0("selection_mode\t", mode_sel),
              paste0("n_selected_modules\t", length(sel)), paste0("collections\t", paste(collections, collapse = ",")),
              paste0(names(opt), "\t", vapply(opt, function(x) if (is.null(x)) "" else paste(x, collapse = ","), "")),
              if (!is.null(gene_sets)) paste0("gene_sets_", names(gene_sets$meta), "\t",
                                              vapply(gene_sets$meta, function(x) paste(x, collapse = ","), "")) else character(0),
              if (!is.null(extended_tss)) paste0("great_", names(attr(extended_tss, "great_params")), "\t",
                                                 unlist(attr(extended_tss, "great_params"))) else character(0),
              if (opt$do_enrichr) paste0("enrichr_status\t", if (isTRUE(enrichr_state$ok)) "ok" else enrichr_state$reason) else character(0))
  write_lines_safe(params, file.path(out_dir, "run_parameters.txt"))
  write_session_info(file.path(out_dir, "session_info.txt"),
                     extra_files = list(gene_set_cache = opt$gene_set_cache, stats_file = stats_file,
                                        module_list_file = module_list_file))
  L("Finished enrichment for ", variant_label)
  invisible(out_dir)
}

# ============================================================
# 10) Enrichr preflight (once), run configs x variants
# ============================================================
enrichr_state <- NULL
if (opt$do_enrichr) {
  enrichr_state <- enrichr_preflight(trimws(strsplit(opt$enrichr_dbs, ",")[[1]]), opt$use_enrichr_background)
  if (!enrichr_state$ok) {
    message("\n*** ENRICHR WILL NOT RUN: ", enrichr_state$reason, " ***\n")
  } else if (length(enrichr_state$dropped) > 0) {
    message("Enrichr: dropping unknown libraries: ", paste(enrichr_state$dropped, collapse = ", "))
  }
}

module_list_files <- resolve_list_input(opt$module_list_file)
trait_list_files <- resolve_list_input(opt$trait_list_file)
configs <- list()
if (!is.null(module_list_files)) {
  for (i in seq_along(module_list_files)) configs[[i]] <- list(stem = names(module_list_files)[i],
                                                             module_list_file = unname(module_list_files[i]), trait_list = NULL)
} else if (!is.null(trait_list_files)) {
  validate_file_exists(opt$stats_file_v1, "--stats_file_v1")
  for (i in seq_along(trait_list_files)) configs[[i]] <- list(stem = names(trait_list_files)[i], module_list_file = NULL,
                                                            trait_list = read_clean_txt_list(unname(trait_list_files[i])))
} else {
  configs[[1]] <- list(stem = "", module_list_file = NULL, trait_list = NULL)
}

# Variant label = the 12a folder name (e.g. v2_exclude_protected_pcs), so logs
# are right even when a v2 folder is passed as --annotation_dir_v1.
variants <- list(v1 = list(dir = opt$annotation_dir_v1, stats = opt$stats_file_v1),
                 v2 = list(dir = opt$annotation_dir_v2, stats = opt$stats_file_v2),
                 v3 = list(dir = opt$annotation_dir_v3, stats = opt$stats_file_v3))
for (k in names(variants)) if (!is.null(variants[[k]]$dir)) variants[[k]]$label <- basename(normalizePath(variants[[k]]$dir))
for (cfg in configs) {
  for (v in variants) {
    if (is.null(v$dir)) next
    if (!is.null(cfg$trait_list)) validate_file_exists(v$stats, paste0("stats file for ", v$label))
    run_enrichment_for_variant(v$dir, v$label, stats_file = v$stats, module_list_file = cfg$module_list_file,
                               trait_list = cfg$trait_list, list_stem = cfg$stem, enrichr_state = enrichr_state)
  }
}

message("Script 12b complete: Enrichment finished")
