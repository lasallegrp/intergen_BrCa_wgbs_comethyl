#!/usr/bin/env Rscript
# =============================================================================
# blood_cell_composition.R
#
# Reference-based blood cell-type deconvolution from Bismark CpG reports
# (coverage2cytosine *.CpG_report.txt.gz) -- no BAMs required.
#
# Steps
#   1. Build windows (+/- --window bp) around 450K/EPIC reference CpGs
#      (sesameData manifests, or a user-supplied TSV for offline nodes).
#   2. For each CpG report: collapse strands, pool M/U counts per window.
#   3. Beta = M / (M+U); windows with coverage < --min_cov set to NA.
#   4. Deconvolve with
#        - EpiDISH RPC  + centDHSbloodDMC.m (B, NK, CD4T, CD8T, Mono, Neutro, Eosino)
#        - Houseman CP  + IDOL-optimized EPIC reference (CD8T, CD4T, NK, Bcell, Mono, Neu)
#   5. Write proportions, per-sample probe-recovery QC, and a method-agreement table/plot.
#
# Usage (example: chr22 test set)
#   pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/tests/scripts/cell_type_decon/blood_cell_composition.R \
#     --input_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/tests/data/processed/cytosine_reports \
#     --genome hg38 --chroms chr22 \
#     --window 100 --min_cov 5 --cores 8 \
#     --outdir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/tests/results/cell_composition_chr22
#  --max_na_frac 0.1
#
# Bismark CpG report columns (no header):
#   chr  pos  strand  count_meth  count_unmeth  context  trinucleotide
# Each CpG appears twice (+ strand C at pos, - strand C at pos+1); the
# '-' strand row is shifted to pos-1 and summed with '+'.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(BiocParallel)
})

opt_list <- list(
  make_option("--input_dir",   type = "character", default = NULL,
              help = "Directory containing CpG reports"),
  make_option("--file_list",   type = "character", default = NULL,
              help = "Optional text file: one CpG-report filename/path per line (overrides --pattern)"),
  make_option("--pattern",     type = "character", default = "CpG_report\\.txt\\.gz$",
              help = "Regex for CpG report files in --input_dir [default %default]"),
  make_option("--genome",      type = "character", default = "hg38",
              help = "hg38 or hg19 -- must match the alignment [default %default]"),
  make_option("--manifest_tsv",type = "character", default = NULL,
              help = "Offline alternative to sesameData: TSV with columns probe, chr, pos"),
  make_option("--chroms",      type = "character", default = NULL,
              help = "Comma-separated chromosomes to keep, e.g. chr22 [default all]"),
  make_option("--window",      type = "integer",   default = 75L,
              help = "Pool CpGs within +/- this many bp of each probe; 0 = probe CpG only [default %default]"),
  make_option("--min_cov",     type = "integer",   default = 5L,
              help = "Minimum pooled reads per window to call a beta [default %default]"),
  make_option("--max_na_frac", type = "double",    default = 0.2,
              help = "Drop reference CpGs missing in more than this fraction of samples [default %default]"),
  make_option("--cores",       type = "integer",   default = 4L,
              help = "Parallel workers for reading reports [default %default]"),
  make_option("--outdir",      type = "character", default = "cell_composition",
              help = "Output directory [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
add_chr <- function(x) ifelse(startsWith(x, "chr"), x, paste0("chr", x))
keep_chroms <- if (is.null(opt$chroms)) NULL else add_chr(trimws(strsplit(opt$chroms, ",")[[1]]))

# -----------------------------------------------------------------------------
# 0. Input files and sample IDs
# -----------------------------------------------------------------------------
if (!is.null(opt$file_list)) {
  files <- trimws(readLines(opt$file_list))
  files <- files[nzchar(files)]
  if (!is.null(opt$input_dir)) {
    rel <- !startsWith(files, "/")
    files[rel] <- file.path(opt$input_dir, files[rel])
  }
} else {
  stopifnot(!is.null(opt$input_dir))
  files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
}
missing <- files[!file.exists(files)]
if (length(missing)) stop("Missing files:\n", paste(missing, collapse = "\n"))
if (!length(files)) stop("No CpG reports found.")

# "4001_merged_name_sorted.deduplicated...CpG_report.txt.gz" -> "4001"
sample_ids <- sub("_merged_name_sorted.*$", "", basename(files))
sample_ids <- sub("\\.CpG_report\\.txt\\.gz$", "", sample_ids)
if (anyDuplicated(sample_ids)) stop("Duplicated sample IDs after parsing filenames.")
names(files) <- sample_ids
msg("Found ", length(files), " CpG reports.")

# -----------------------------------------------------------------------------
# 1. References and probe coordinates
# -----------------------------------------------------------------------------
data("centDHSbloodDMC.m", package = "EpiDISH", envir = environment())
ref_epidish <- centDHSbloodDMC.m

data("IDOLOptimizedCpGs.compTable", package = "FlowSorted.Blood.EPIC", envir = environment())
idol_cells <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
if (!all(idol_cells %in% colnames(IDOLOptimizedCpGs.compTable)))
  stop("Unexpected IDOL compTable columns: ",
       paste(colnames(IDOLOptimizedCpGs.compTable), collapse = ", "))
ref_idol <- as.matrix(IDOLOptimizedCpGs.compTable[, idol_cells])

ref_probes <- union(rownames(ref_epidish), rownames(ref_idol))
msg("Reference CpGs: EpiDISH=", nrow(ref_epidish), ", IDOL=", nrow(ref_idol),
    ", union=", length(ref_probes))

if (!is.null(opt$manifest_tsv)) {
  man <- fread(opt$manifest_tsv)
  stopifnot(all(c("probe", "chr", "pos") %in% names(man)))
  man <- man[, .(probe, chr = add_chr(as.character(chr)), pos = as.integer(pos))]
} else {
  # sesame manifests: CpG coordinate per probe; EPIC first, HM450 fills
  # 450K-only probes used by centDHSbloodDMC.m. 1-bp convention differences
  # are absorbed by the window.
  man <- rbindlist(lapply(c("EPIC", "HM450"), function(p) {
    gr <- sesameData::sesameData_getManifestGRanges(p, genome = opt$genome)
    data.table(probe = names(gr), chr = as.character(GenomicRanges::seqnames(gr)),
               pos = GenomicRanges::start(gr))
  }))
  man <- unique(man, by = "probe")
}
man <- man[probe %in% ref_probes & !is.na(pos)]
if (!is.null(keep_chroms)) man <- man[chr %in% keep_chroms]
msg("Reference CpGs with coordinates", if (!is.null(keep_chroms)) " on selected chromosomes", ": ",
    nrow(man), " (EpiDISH ", sum(man$probe %in% rownames(ref_epidish)),
    ", IDOL ", sum(man$probe %in% rownames(ref_idol)), ")")
if (nrow(man) == 0) stop("No reference CpGs on the selected chromosomes / genome build.")

win <- man[, .(probe, chr, start = pos - opt$window, end = pos + 1L + opt$window)]
setkey(win, chr, start, end)

# -----------------------------------------------------------------------------
# 2. Pool counts per window, per sample
# -----------------------------------------------------------------------------
summarize_report <- function(f, win, keep_chroms) {
  library(data.table)
  dt <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = 1:5,
              col.names = c("chr", "pos", "strand", "M", "U"),
              colClasses = list(character = c(1, 3), integer = c(2, 4, 5)))
  dt <- dt[M + U > 0L]
  if (nrow(dt) && !startsWith(dt$chr[1], "chr")) dt[, chr := paste0("chr", chr)]
  if (!is.null(keep_chroms)) dt <- dt[chr %in% keep_chroms]
  dt[strand == "-", pos := pos - 1L]
  dt <- dt[, .(M = sum(M), U = sum(U)), by = .(chr, pos)]

  qc <- data.table(n_cpg_covered = nrow(dt),
                   mean_depth    = if (nrow(dt)) mean(dt$M + dt$U) else NA_real_,
                   median_depth  = if (nrow(dt)) as.numeric(median(dt$M + dt$U)) else NA_real_)

  j <- win[dt, on = .(chr, start <= pos, end >= pos), nomatch = NULL,
           allow.cartesian = TRUE, .(probe = x.probe, M = i.M, U = i.U)]
  j <- j[, .(M = sum(M), U = sum(U), n_cpg = .N), by = probe]
  list(counts = j, qc = qc)
}

bp <- if (opt$cores > 1) MulticoreParam(workers = opt$cores) else SerialParam()
msg("Reading reports with ", opt$cores, " worker(s)...")
res <- bplapply(files, summarize_report, win = win, keep_chroms = keep_chroms, BPPARAM = bp)

probes <- man$probe
Mmat <- Cmat <- matrix(0L, nrow = length(probes), ncol = length(files),
                       dimnames = list(probes, sample_ids))
for (s in sample_ids) {
  cnt <- res[[s]]$counts
  idx <- match(cnt$probe, probes)
  Mmat[idx, s] <- cnt$M
  Cmat[idx, s] <- cnt$M + cnt$U
}
beta <- Mmat / Cmat
beta[Cmat < opt$min_cov] <- NA

depth_qc <- rbindlist(lapply(sample_ids, function(s) cbind(sample = s, res[[s]]$qc)))

fwrite(data.table(probe = probes, beta), file.path(opt$outdir, "reference_cpg_betas.tsv.gz"), sep = "\t")
fwrite(data.table(probe = probes, Cmat), file.path(opt$outdir, "reference_cpg_coverage.tsv.gz"), sep = "\t")

# -----------------------------------------------------------------------------
# 3. Per-sample probe recovery QC
# -----------------------------------------------------------------------------
recov <- function(ref) {
  p <- intersect(rownames(ref), probes)
  if (!length(p)) return(rep(0L, length(sample_ids)))
  colSums(!is.na(beta[p, , drop = FALSE]))
}
qc <- depth_qc[, `:=`(
  epidish_ref_total     = nrow(ref_epidish),
  epidish_ref_recovered = recov(ref_epidish),
  idol_ref_total        = nrow(ref_idol),
  idol_ref_recovered    = recov(ref_idol)
)]
fwrite(qc, file.path(opt$outdir, "qc_probe_recovery.tsv"), sep = "\t")
msg("Median recovered reference CpGs/sample: EpiDISH ", median(qc$epidish_ref_recovered),
    "/", nrow(ref_epidish), ", IDOL ", median(qc$idol_ref_recovered), "/", nrow(ref_idol))

# -----------------------------------------------------------------------------
# 4. Deconvolution
# -----------------------------------------------------------------------------
prep_matrix <- function(ref, label) {
  p <- intersect(rownames(ref), probes)
  b <- beta[p, , drop = FALSE]
  keep <- rowMeans(is.na(b)) <= opt$max_na_frac
  b <- b[keep, , drop = FALSE]
  if (nrow(b)) {
    med <- apply(b, 1, median, na.rm = TRUE)
    na_idx <- which(is.na(b), arr.ind = TRUE)
    b[na_idx] <- med[na_idx[, 1]]
  }
  n_min <- 3 * ncol(ref)
  msg(label, ": ", nrow(b), " reference CpGs usable after NA filter")
  if (nrow(b) < n_min) {
    warning(label, ": only ", nrow(b), " usable CpGs (< ", n_min,
            "); skipping. Expected for single-chromosome test data.", call. = FALSE)
    return(NULL)
  }
  if (nrow(b) < 0.5 * nrow(ref))
    warning(label, ": using ", nrow(b), "/", nrow(ref),
            " reference CpGs; treat estimates as provisional.", call. = FALSE)
  list(beta = b, ref = ref[rownames(b), , drop = FALSE])
}

results <- list()

d <- prep_matrix(ref_epidish, "EpiDISH RPC")
if (!is.null(d)) {
  est <- EpiDISH::epidish(beta.m = d$beta, ref.m = d$ref, method = "RPC")$estF
  results$epidish <- data.table(sample = rownames(est), est, n_cpg_used = nrow(d$beta))
  fwrite(results$epidish, file.path(opt$outdir, "cellprop_epidish_RPC.tsv"), sep = "\t")
}

d <- prep_matrix(ref_idol, "IDOL Houseman CP")
if (!is.null(d)) {
  est <- FlowSorted.Blood.EPIC::projectCellType_CP(d$beta, d$ref,
                                                   nonnegative = TRUE, lessThanOne = FALSE)
  # CP estimates are not constrained to sum to 1; keep raw and add normalized.
  est_norm <- est / rowSums(est)
  results$idol <- data.table(sample = rownames(est), est,
                             sum_raw = rowSums(est), n_cpg_used = nrow(d$beta))
  fwrite(results$idol, file.path(opt$outdir, "cellprop_idol_CP_raw.tsv"), sep = "\t")
  fwrite(data.table(sample = rownames(est_norm), est_norm),
         file.path(opt$outdir, "cellprop_idol_CP_normalized.tsv"), sep = "\t")
  results$idol_norm <- data.table(sample = rownames(est_norm), est_norm)
}

# -----------------------------------------------------------------------------
# 5. Method agreement (shared cell types)
# -----------------------------------------------------------------------------
if (!is.null(results$epidish) && !is.null(results$idol_norm)) {
  map <- c(B = "Bcell", NK = "NK", CD4T = "CD4T", CD8T = "CD8T", Mono = "Mono", Neutro = "Neu")
  long_e <- melt(results$epidish[, c("sample", names(map)), with = FALSE],
                 id.vars = "sample", variable.name = "cell", value.name = "epidish")
  long_e[, cell := map[as.character(cell)]]
  long_i <- melt(results$idol_norm[, c("sample", unname(map)), with = FALSE],
                 id.vars = "sample", variable.name = "cell", value.name = "idol")
  long_i[, cell := as.character(cell)]
  cmp <- merge(long_e, long_i, by = c("sample", "cell"))
  agree <- cmp[, .(pearson_r = cor(epidish, idol), spearman_rho = cor(epidish, idol, method = "spearman"),
                   mean_abs_diff = mean(abs(epidish - idol)),
                   mean_epidish = mean(epidish), mean_idol = mean(idol)), by = cell]
  fwrite(agree, file.path(opt$outdir, "method_agreement.tsv"), sep = "\t")
  print(agree)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    g <- ggplot(cmp, aes(idol, epidish)) +
      geom_abline(linetype = 2, colour = "grey60") +
      geom_point(size = 1, alpha = 0.6) +
      facet_wrap(~cell, scales = "free") +
      labs(x = "IDOL / Houseman CP (normalized)", y = "EpiDISH RPC",
           title = "Blood cell-type estimates: method agreement") +
      theme_bw(base_size = 10)
    ggsave(file.path(opt$outdir, "method_agreement.pdf"), g, width = 8, height = 5.5)
  }
}

writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))
msg("Done. Outputs in ", normalizePath(opt$outdir))