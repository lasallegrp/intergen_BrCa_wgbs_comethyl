# Blood cell-type deconvolution from WGBS (`blood_cell_composition.R`)

Estimates blood cell-type proportions from Bismark CpG reports
(`*.CpG_report.txt.gz`). BAM files are not needed.

## Method

1. **Reference CpGs.** These are the CpGs in two established blood references:
   - **EpiDISH** `centDHSbloodDMC.m`: 333 CpGs; B, NK, CD4T, CD8T, Mono, Neutro, Eosino
   - **IDOL** `IDOLOptimizedCpGs.compTable` (FlowSorted.Blood.EPIC): 450 CpGs; CD8T, CD4T, NK, Bcell, Mono, Neu

   Probe coordinates come from the sesameData EPIC and HM450 manifests for the
   chosen genome build, or from a TSV given with `--manifest_tsv`.
2. **Strand collapse.** Each CpG appears twice in a CpG report. The `-` strand
   row is shifted to `pos - 1` and its counts are added to the `+` strand row.
3. **Window pooling.** Methylated and unmethylated counts from all CpGs within
   ±`--window` bp of each reference CpG are summed. At ~8× depth a single CpG
   is too noisy (SD ≈ 0.18 at β = 0.5). Pooling 3–10 neighbouring CpGs, which
   share methylation state in blood, brings this down to about 0.06–0.10.
4. **Beta calling.** β = M / (M + U). Windows with fewer than `--min_cov`
   pooled reads are set to NA. CpGs missing in more than `--max_na_frac` of
   samples are dropped; the remaining NAs are filled with the CpG's median
   across samples.
5. **Deconvolution.**
   - EpiDISH robust partial correlations (RPC)
   - Houseman constrained projection (CP) with the IDOL reference. The raw
     estimates are not forced to sum to 1, so both raw and row-normalized
     values are written.

   A method is skipped with a warning if fewer than 3 × (number of cell types)
   reference CpGs are usable. This is expected for single-chromosome test data.
6. **Agreement.** Estimates from the two methods are compared per shared cell
   type: Pearson, Spearman and mean absolute difference, plus a scatter plot.

## Usage

```bash
pixi run Rscript blood_cell_composition.R \
  --input_dir /path/to/cytosine_reports \
  --genome hg38 \
  --window 100 --min_cov 15 --max_na_frac 0.1 \
  --cores 8 \
  --outdir results/cell_composition
```

| Option | Default | Notes |
|---|---|---|
| `--input_dir` | – | Directory of CpG reports |
| `--file_list` | – | Optional list of report files (overrides `--pattern`) |
| `--pattern` | `CpG_report\.txt\.gz$` | Regex for finding reports |
| `--genome` | `hg38` | Must match the alignment build |
| `--manifest_tsv` | – | Offline probe coordinates (`probe`, `chr`, `pos`) |
| `--chroms` | all | e.g. `chr22` for testing |
| `--window` | 75 | ± bp pooled around each probe; 0 = probe CpG only |
| `--min_cov` | 5 | Minimum pooled reads per window |
| `--max_na_frac` | 0.2 | Maximum fraction of samples allowed to be missing per CpG |
| `--cores` | 4 | Parallel workers for reading reports |
| `--outdir` | `cell_composition` | Output directory |

Sample IDs are taken from file names: the part before `_merged_name_sorted`
(for example `4001`).

## Outputs

| File | Contents |
|---|---|
| `qc_probe_recovery.tsv` | Per sample: CpGs covered, mean and median depth, reference CpGs recovered for each method |
| `reference_cpg_betas.tsv.gz` / `reference_cpg_coverage.tsv.gz` | Pooled β values and coverage per reference CpG |
| `cellprop_epidish_RPC.tsv` | EpiDISH proportions |
| `cellprop_idol_CP_raw.tsv` / `cellprop_idol_CP_normalized.tsv` | IDOL/Houseman proportions |
| `method_agreement.tsv` / `.pdf` | Agreement between methods |
| `sessionInfo.txt` | Package versions |

## Setup

pixi dependencies (channels: conda-forge, bioconda):
`r-base`, `r-optparse`, `r-data.table`, `r-ggplot2`,
`bioconductor-biocparallel`, `bioconductor-genomicranges`,
`bioconductor-epidish`, `bioconductor-sesamedata`,
`bioconductor-flowsorted.blood.epic`.

Download the sesame data once, on a node with internet access. The cache
location is set in `~/.bashrc` and reused by every project:

```bash
export EXPERIMENT_HUB_CACHE=/quobyte/lasallegrp/home/George/.cache/ExperimentHub
pixi run Rscript -e 'sesameData::sesameDataCache()'   # all sesameData files, a few GB
# check: both should print probe counts (~865k EPIC, ~485k HM450)
pixi run Rscript -e 'library(sesameData); for (p in c("EPIC","HM450")) print(length(sesameData_getManifestGRanges(p, genome = "hg38")))'
```

The sesameData version in this environment (R 4.3) does not accept
`sesameDataCache("EPIC.address")`. To download only the two required files,
fetch them through ExperimentHub instead (`query(ExperimentHub(), c("sesameData", "EPIC.address"))`),
or skip sesame entirely with `--manifest_tsv`.

## Notes and caveats

- **CHDS data:** 173 samples, ~8× mean CpG depth (range 6.2–11.3×), 77% of
  CpGs at ≥5×. Recommended settings: `--window 100 --min_cov 15`. Relax to
  `--min_cov 10` if the lowest-depth samples (4148, 4042, 4102) lose more than
  10% of reference CpGs.
- **chr22 test:** this only checks that the pipeline runs end to end. Only a
  handful of reference CpGs lie on chr22, so no proportions are estimated.
- **Depth check:** correlate each cell proportion with sample mean depth. If
  they are correlated, add depth as a covariate in downstream models.
- **Collinearity:** proportions sum to ~1. When adjusting for them, drop one
  cell type or use the first 1–2 PCs of the proportions.
- **Interpretation:** these are reference-based estimates from probe-mapped
  WGBS, not measured cell counts.