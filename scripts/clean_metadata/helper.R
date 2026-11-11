# ================================================================
# helper.R
#
# Shared helper functions for the reproducible comethyl pipeline.
#
# FUNCTIONS INCLUDED
#   - TRAIT / SAMPLE INFO (scripts 04-09):
#       readTraitFile(), readSampleInfo(), resolveTraits(),
#       collect_set_files(), read_trait_set_file(), get_set_name()
#   - PC-TRAIT AND ME-TRAIT PLOTS / OUTPUTS:
#       .standardize_pc_trait_table(), plotPCTrait(),
#       plotTraitDendrogramFromPC(), plotMEtraitCor(),
#       save_me_trait_method_outputs(), validate_modules_object()
#   - ME-TRAIT PAIR PLOTS:
#       auto_ylim(), is_binary_like(), relabel_binary_to_yesno(),
#       make_me_trait_subtitle(), add_plot_subtitle(), plotMEtraitViolin()
#   - CLI / LOGGING / VALIDATION / REPRODUCIBILITY:
#       get_arg(), trim_or_null(), parse_bool(), safe_dir_create(),
#       timestamp_now(), append_log(), write_lines_safe(),
#       write_log_lines(), write_vector_file(), stop_if_missing(),
#       validate_file_exists(), validate_regions_df(),
#       extract_regions_from_module_object(), write_session_info(),
#       file_checksum(), resolve_annotation_packages()
#   - GREAT REGION -> GENE (DEFAULT gene assignment; matches
#     comethyl::annotateModule(); offline via rGREAT; used by 12a and 12b):
#       great_parse_tss_source(), great_build_extended_tss(),
#       great_region_gene_associations(), great_collapse_associations(),
#       gene_info_from_orgdb()
#   - ANNOTATR GENIC / ENHANCER / CpG CONTEXT (matches
#     comethyl::annotateModule(); used by 12a):
#       annotatr_category_table(), annotatr_default_priority(),
#       annotatr_genome_seqinfo(), build_annotatr_annotations()
#       [+ build_cpg_annotations_annotatr_local(),
#        build_enhancer_annotations_local()], annotate_regions_annotatr(),
#       moduleAnnotatrEnrichment(), summarizeAnnotationDistribution(),
#       plotAnnotationDistribution(), plotAnnotatrEnrichmentHeatmap(),
#       plotAnnotatrEnrichmentForest() [+ annot_* plot utilities]
#   - GENE-SET CACHE + FAST GREAT TEST (12b):
#       parse_collections(), build_gene_set_cache(), load_gene_set_cache(),
#       symbols_to_entrez(), great_fast_prepare(), great_fast_test()
#   - NEAREST-GENE PATH (12a --gene_assignment nearest; NOT comethyl's rule):
#       offline_nearest_gene()
#   - LEGACY ChIPseeker/DMRichR-style annotation + enrichment (pre-annotatr
#     12a; kept so earlier outputs stay reproducible; NOT called by the
#     current 12a/12b):
#       annotate_offline_only(), annotate_regions_safe(),
#       annotate_genomic_location(), build_cpg_annotations(),
#       classify_cpg_flags(), derive_cpg_context_label(),
#       moduleGenicEnrichment(), moduleCpGEnrichment(),
#       plotModuleAnnotationEnrichment()
#
# NOTES
#   - These helpers are generic and do not include project-specific
#     trait recoding rules.
#   - sample_info is assumed to be analysis-ready before entering the
#     pipeline.
#   - Annotation helpers call GenomicRanges/IRanges/AnnotationDbi/
#     GenomeInfoDb/S4Vectors/comethyl functions via pkg::fn()
#     namespacing, so those packages just need to be installed
#     (not attached) in the environment sourcing this file.
#
# CHANGE LOG
#   - Merged helper (10) [annotatr / offline-GREAT / gene-set cache
#     stack used by 12a/12b] with helper (11) [restored 7b5ab0b
#     trait-handling fixes, lost in fac40d7]. From helper (11):
#       * readSampleInfo(): readxl for .xlsx, .name_repair = "minimal",
#         so names stay exactly as in Excel ("F1 Breast Cancer", not
#         "F1.Breast.Cancer"); date columns come back as dates and are
#         dropped as non-numeric. First Excel column = sample IDs.
#       * resolveTraits(): case/space/punctuation-insensitive matching;
#         returns the ACTUAL sample_info column name; stops on
#         ambiguous collisions.
#       * plotMEtraitCor(): publication defaults (18/18/18, color-bar
#         label size 4); x-axis hidden whenever showColorBar = TRUE.
#     All annotation/enrichment/GREAT/gene-set code = helper (10).
#   - Documentation pass: every function now has a header comment.
#     Corrected stale comments: GREAT (not nearest gene) is the default
#     gene assignment in 12a, matching comethyl; the annotate_offline_only/
#     annotate_regions_safe path is legacy; plotPCTrait() describes its
#     real arguments (subsetTraits / topCells), and its messages are
#     tagged [plotPCTrait] instead of [plotPCTrait].
# ================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(readxl)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(scales)
})

# ------------------------------------------------------------
# readTraitFile
# ------------------------------------------------------------
# Reads a plain text file with one trait name per line.
# Removes blank lines, trims whitespace, and returns unique names.
#
# Example file:
#   BrCAF1
#   PCB138
#   DDE
# ------------------------------------------------------------
readTraitFile <- function(file, verbose = TRUE) {
  if (is.null(file)) return(NULL)

  if (!file.exists(file)) {
    stop("Trait file not found: ", file)
  }

  traits <- readLines(file, warn = FALSE)
  traits <- trimws(traits)
  traits <- traits[traits != ""]
  traits <- unique(traits)

  if (verbose) {
    message("[readTraitFile] Loaded ", length(traits), " traits from: ", file)
  }

  return(traits)
}

# ------------------------------------------------------------
# readSampleInfo
# ------------------------------------------------------------
# Reads sample information from .xlsx, .csv, or .tsv.
#
# Behavior:
#   - .xlsx/.xls: read with readxl (sheet 1). The FIRST column is
#     taken as sample IDs and becomes the rownames; it is removed
#     from the returned columns. Column names are kept exactly as
#     written in Excel (.name_repair = "minimal").
#   - If sample_id_col is provided, that column is used as rownames.
#   - If sample_id_col is not provided, existing rownames are used.
#   - If rownames fallback is used, the function prints the first few
#     rownames so the user can confirm they are correct.
#
# Returns:
#   data.frame with sample IDs in rownames
# ------------------------------------------------------------
readSampleInfo <- function(file,
                           sample_id_col = NULL,
                           verbose = TRUE) {
  if (!file.exists(file)) {
    stop("Sample info file not found: ", file)
  }

  ext <- tolower(tools::file_ext(file))

  if (verbose) {
    message("[readSampleInfo] Reading sample info: ", file)
  }

  if (ext %in% c("xlsx", "xls")) {
    raw_df <- readxl::read_excel(
      path = file,
      sheet = 1,
      .name_repair = "minimal"
    )

    raw_df <- as.data.frame(
      raw_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # First Excel column contains the sample IDs / original rownames
    sample_ids <- trimws(as.character(raw_df[[1]]))

    if (any(is.na(sample_ids) | sample_ids == "")) {
      stop("First Excel column contains missing or blank sample IDs.")
    }

    if (anyDuplicated(sample_ids)) {
      dup_ids <- unique(sample_ids[duplicated(sample_ids)])
      stop(
        "Duplicate sample IDs found in first Excel column. Examples: ",
        paste(head(dup_ids, 10), collapse = ", ")
      )
    }

    df <- raw_df[, -1, drop = FALSE]
    rownames(df) <- sample_ids

  } else if (ext == "csv") {
    df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)

  } else if (ext %in% c("tsv", "txt")) {
    df <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)

  } else {
    stop("Unsupported sample info format: .", ext,
         " (supported: .xlsx, .csv, .tsv, .txt)")
  }

  if (nrow(df) == 0) {
    stop("Sample info file is empty: ", file)
  }

  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(df)) {
      stop("sample_id_col not found in sample info: ", sample_id_col)
    }

    sample_ids <- as.character(df[[sample_id_col]])
    sample_ids <- trimws(sample_ids)

    if (any(is.na(sample_ids) | sample_ids == "")) {
      stop("sample_id_col contains missing or blank sample IDs: ", sample_id_col)
    }

    if (anyDuplicated(sample_ids)) {
      dup_ids <- unique(sample_ids[duplicated(sample_ids)])
      stop("Duplicate sample IDs found in sample_id_col '", sample_id_col,
           "'. Example duplicates: ", paste(head(dup_ids, 10), collapse = ", "))
    }

    rownames(df) <- sample_ids

    if (verbose) {
      message("[readSampleInfo] Using sample ID column: ", sample_id_col)
      message("[readSampleInfo] First 6 sample IDs:")
      print(utils::head(rownames(df)))
    }

  } else {
    rn <- rownames(df)

    if (is.null(rn) || length(rn) == 0 || any(rn == "")) {
      stop(
        "No sample_id_col provided and rownames are not available/usable. ",
        "Please provide --sample_id_col."
      )
    }

    if (anyDuplicated(rn)) {
      dup_ids <- unique(rn[duplicated(rn)])
      stop("Duplicate sample IDs found in rownames. Example duplicates: ",
           paste(head(dup_ids, 10), collapse = ", "))
    }

    if (verbose) {
      message("[readSampleInfo] No sample_id_col provided; using existing rownames.")
      message("[readSampleInfo] First 6 rownames used as sample IDs:")
      print(utils::head(rn))
    }
  }

  return(df)
}

# ------------------------------------------------------------
# resolveTraits
# ------------------------------------------------------------
# Compares requested trait names against available trait columns.
#
# Matching ignores case, whitespace, and punctuation, so a trait file
# entry "PCB 138", "pcb.138" or "PCB-138" matches a column "PCB138".
# The ACTUAL column name from available_traits is returned in $found.
# Stops if two available columns collapse to the same key (ambiguous).
#
# Returns a list with:
#   requested
#   found
#   missing
# ------------------------------------------------------------
resolveTraits <- function(requested_traits,
                          available_traits,
                          label = "traits",
                          verbose = TRUE) {

  if (is.null(requested_traits)) {
    return(list(
      requested = character(0),
      found = character(0),
      missing = character(0)
    ))
  }

  requested_traits <- unique(trimws(as.character(requested_traits)))
  requested_traits <- requested_traits[nzchar(requested_traits)]

  available_traits <- as.character(available_traits)

  # Comparison key that ignores differences such as:
  # spaces vs dots, parentheses, colon, slash, hyphen, case, etc.
  trait_key <- function(x) {
    x <- trimws(tolower(as.character(x)))
    gsub("[[:space:][:punct:]]+", "", x)
  }

  requested_keys <- trait_key(requested_traits)
  available_keys <- trait_key(available_traits)

  # Stop if two available columns collapse to the same matching key.
  duplicated_available <- unique(available_keys[duplicated(available_keys)])
  if (length(duplicated_available) > 0) {
    ambiguous <- lapply(duplicated_available, function(k) {
      available_traits[available_keys == k]
    })

    stop(
      "Ambiguous trait-name matching after punctuation normalization. ",
      "Examples:\n",
      paste(
        vapply(
          ambiguous,
          function(x) paste("  ", paste(x, collapse = " | ")),
          character(1)
        ),
        collapse = "\n"
      )
    )
  }

  found <- character(0)
  missing <- character(0)

  for (i in seq_along(requested_traits)) {
    hit <- which(available_keys == requested_keys[i])

    if (length(hit) == 1) {
      # Return the ACTUAL name in the stats/sample-info table.
      found <- c(found, available_traits[hit])
    } else {
      missing <- c(missing, requested_traits[i])
    }
  }

  found <- unique(found)
  missing <- unique(missing)

  if (verbose) {
    message(
      "[resolveTraits] ", label,
      ": requested = ", length(requested_traits),
      ", found = ", length(found),
      ", missing = ", length(missing)
    )

    if (length(missing) > 0) {
      message(
        "[resolveTraits] Missing ", label, ": ",
        paste(head(missing, 10), collapse = ", ")
      )
    }
  }

  list(
    requested = requested_traits,
    found = found,
    missing = missing
  )
}

# ------------------------------------------------------------
# .standardize_pc_trait_table
# ------------------------------------------------------------
# Internal helper to standardize PC/trait stats table.
# Requires:
#   - PC column: "module" or "PC"
#   - trait column: "trait"
#   - p column: "p"
#   - effect size column: "cor" or "bicor"
# ------------------------------------------------------------
.standardize_pc_trait_table <- function(pc_trait_stats,
                                        cor_column = c("bicor", "cor"),
                                        p_column = "p") {
  cor_column <- match.arg(cor_column)

  df <- as.data.frame(pc_trait_stats, stringsAsFactors = FALSE)

  pc_col <- NULL
  if ("module" %in% colnames(df)) {
    pc_col <- "module"
  } else if ("PC" %in% colnames(df)) {
    pc_col <- "PC"
  } else {
    stop("pc_trait_stats must contain either 'module' or 'PC' column.")
  }

  if (!"trait" %in% colnames(df)) {
    stop("pc_trait_stats must contain a 'trait' column.")
  }

  if (!p_column %in% colnames(df)) {
    stop("pc_trait_stats must contain p-value column: ", p_column)
  }

  if (!cor_column %in% colnames(df)) {
    stop("pc_trait_stats must contain effect-size column: ", cor_column)
  }

  out <- df %>%
    dplyr::transmute(
      PC = as.character(.data[[pc_col]]),
      trait = as.character(.data[["trait"]]),
      effect = as.numeric(.data[[cor_column]]),
      p = as.numeric(.data[[p_column]])
    )

  out <- out %>%
    dplyr::filter(!is.na(PC), !is.na(trait), !is.na(effect), !is.na(p))

  if (nrow(out) == 0) {
    stop("No valid rows remain in pc_trait_stats after standardization.")
  }

  return(out)
}

# ------------------------------------------------------------
# plotPCTrait
# ------------------------------------------------------------
# Heatmap of PC-trait (or ME-trait) association statistics
# (fill = bicor/cor, star or p-value label where p < p_cut), with an
# optional module colour strip when column names are R colours.
#
# Input: data.frame with module, trait, p and bicor (or cor) columns.
# Row/column subsetting:
#   subsetTraits = character vector -> plot only these traits
#   topCells     = N                -> keep the N smallest-p cells
#   (both NULL                      -> plot everything)
# Axis labels are hidden past maxXLabels / maxYLabels unless
# force_x_labels / force_y_labels. autoSize scales the figure to the
# number of rows/columns within [min, max] width/height.
# Returns the (combined) plot; saves to `file` when save = TRUE.
# ------------------------------------------------------------
plotPCTrait <- function(MEtraitCor,
                               moduleOrder = 1:length(unique(MEtraitCor$module)),
                               traitOrder  = 1:length(unique(MEtraitCor$trait)),

                               subsetTraits = NULL,
                               topCells = NULL,
                               p_cut = 0.05,

                               label_mode = c("star", "p", "none"),
                               label.size = 8,
                               label.nudge_y = 0,

                               maxXLabels = 40,
                               maxYLabels = 80,
                               force_x_labels = FALSE,   # NEW
                               force_y_labels = FALSE,   # NEW
                               wrapY = 30,

                               colors = WGCNA::blueWhiteRed(100, gamma = 0.9),
                               limit = NULL,

                               base_size = 11,
                               axis.text.size = 8,
                               legend.text.size = 9,
                               legend.title.size = 10,
                               legend.position = "right",

                               showColorBar = TRUE,
                               showColorBarLabels = TRUE,
                               colorBarLabelPos = c("inside","below"),
                               colorBarLabelSize = 2.8,
                               colorBarLabelAngle = 90,
                               colorBarRelHeight = 0.10,
                               colColorMargins = c(-0.7, 4.21, 1.2, 11.07),
                               syncWidths = TRUE,
                               autoContrastLabels = TRUE,

                               save = TRUE,
                               file = "Heatmap.pdf",
                               autoSize = TRUE,
                               cell_in = 0.16,
                               min_width = 6,
                               min_height = 5,
                               max_width = 11,
                               max_height = 11,
                               width = 8,
                               height = 7,
                               dpi = 600,
                               verbose = TRUE) {

  label_mode <- match.arg(label_mode)
  colorBarLabelPos <- match.arg(colorBarLabelPos)

  stopifnot(all(c("module","trait","p") %in% colnames(MEtraitCor)))

  corType <- if ("bicor" %in% colnames(MEtraitCor)) "bicor" else if ("cor" %in% colnames(MEtraitCor)) "cor" else
    stop("[plotPCTrait] Need 'cor' or 'bicor' column")

  if (verbose) message("[plotPCTrait] corType = ", corType)

  if (!is.null(subsetTraits)) {
    MEtraitCor <- MEtraitCor[MEtraitCor$trait %in% subsetTraits, , drop = FALSE]
  }

  if (!is.null(topCells)) {
    MEtraitCor <- MEtraitCor[order(MEtraitCor$p), , drop = FALSE]
    MEtraitCor <- head(MEtraitCor, topCells)
  }

  # --- safer/stable releveling ---
  mod_levels <- unique(as.character(MEtraitCor$module))
  trt_levels <- unique(as.character(MEtraitCor$trait))

  moduleOrder <- moduleOrder[moduleOrder <= length(mod_levels)]
  traitOrder  <- traitOrder[traitOrder  <= length(trt_levels)]

  MEtraitCor$module <- factor(MEtraitCor$module, levels = mod_levels[moduleOrder])
  MEtraitCor$trait  <- factor(MEtraitCor$trait,  levels = trt_levels[rev(traitOrder)])

  MEtraitCor$significant <- (MEtraitCor$p < p_cut) & !is.na(MEtraitCor$p)

  if (wrapY > 0) {
    MEtraitCor$trait <- factor(
      stringr::str_wrap(as.character(MEtraitCor$trait), width = wrapY),
      levels = stringr::str_wrap(levels(MEtraitCor$trait), width = wrapY)
    )
  }

  corData <- MEtraitCor[[corType]]
  if (is.null(limit)) limit <- max(abs(corData), na.rm = TRUE)
  titleName <- paste0(toupper(substring(corType, 1, 1)), substring(corType, 2))

  nX <- length(levels(MEtraitCor$module))
  nY <- length(levels(MEtraitCor$trait))

  if (autoSize) {
    width  <- max(min_width,  min(max_width,  nX * cell_in + 2.5))
    height <- max(min_height, min(max_height, nY * cell_in + 2.5))
    if (verbose) message(sprintf("[plotPCTrait] autoSize -> width=%.2f, height=%.2f (nX=%d,nY=%d)",
                                 width, height, nX, nY))
  }

  # NEW: override suppression when you *need* labels (PCs)
  show_x_labels <- force_x_labels || (nX <= maxXLabels)
  show_y_labels <- force_y_labels || (nY <= maxYLabels)

  heatmap <- ggplot2::ggplot(MEtraitCor, ggplot2::aes(x = module, y = trait)) +
    ggplot2::geom_tile(ggplot2::aes(fill = corData), linewidth = 0) +
    ggplot2::scale_fill_gradientn(titleName, colors = colors, limits = c(-limit, limit)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = if (show_x_labels) ggplot2::element_text(size = axis.text.size, angle = 90, vjust = 0.5) else ggplot2::element_blank(),
      axis.text.y = if (show_y_labels) ggplot2::element_text(size = axis.text.size) else ggplot2::element_blank(),
      axis.ticks.x = if (show_x_labels) ggplot2::element_line() else ggplot2::element_blank(),
      axis.ticks.y = if (show_y_labels) ggplot2::element_line() else ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.text  = ggplot2::element_text(size = legend.text.size),
      legend.title = ggplot2::element_text(size = legend.title.size)
    )

  if (label_mode != "none") {
    if (label_mode == "p") {
      heatmap <- heatmap +
        ggplot2::geom_text(
          ggplot2::aes(label = ifelse(significant, format(p, digits = 1, scientific = TRUE), "")),
          size = label.size, nudge_y = label.nudge_y
        )
    } else {
      heatmap <- heatmap +
        ggplot2::geom_text(
          ggplot2::aes(label = ifelse(significant, "*", "")),
          size = label.size, nudge_y = label.nudge_y
        )
    }
  }

  # module color strip (unchanged)
  uniqMods   <- levels(MEtraitCor$module)
  colModules <- all(uniqMods %in% grDevices::colors())

  colColors <- NULL
  if (colModules && showColorBar) {
    dfBar <- data.frame(module = factor(uniqMods, levels = uniqMods), y = 1)

    if (showColorBarLabels) {
      if (autoContrastLabels) {
        lum <- function(cols) {
          rgb <- grDevices::col2rgb(cols) / 255
          drop(0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ])
        }
        dfBar$labcol <- ifelse(lum(as.character(dfBar$module)) > 0.5, "black", "white")
      } else dfBar$labcol <- "black"
    }

    baseBar <- ggplot2::ggplot(dfBar, ggplot2::aes(x = module, y = y, fill = module)) +
      ggplot2::geom_tile(width = 0.98, height = 0.8) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none",
                     plot.margin = grid::unit(colColorMargins, "lines"))

    if (showColorBarLabels) {
      if (colorBarLabelPos == "inside") {
        baseBar <- baseBar +
          ggplot2::geom_text(ggplot2::aes(label = as.character(module), colour = labcol),
                             angle = colorBarLabelAngle, vjust = 0.5, size = colorBarLabelSize) +
          ggplot2::scale_colour_identity()
      } else {
        baseBar <- baseBar +
          ggplot2::geom_text(ggplot2::aes(y = 0.05, label = as.character(module), colour = labcol),
                             angle = 0, vjust = 1, size = colorBarLabelSize) +
          ggplot2::scale_colour_identity() +
          ggplot2::ylim(0, 1.35)
      }
    }

    colColors <- baseBar
  }

  combined <- heatmap
  if (!is.null(colColors)) {
    if (syncWidths) {
      g_heat  <- ggplot2::ggplotGrob(heatmap)
      g_strip <- ggplot2::ggplotGrob(colColors)
      g_strip$widths <- g_heat$widths
      combined <- gridExtra::arrangeGrob(
        g_heat, g_strip, ncol = 1,
        heights = grid::unit.c(grid::unit(1, "null"),
                               grid::unit(colorBarRelHeight, "null"))
      )
    } else {
      combined <- cowplot::plot_grid(heatmap, colColors, nrow = 2,
                                     rel_heights = c(1, colorBarRelHeight))
    }
  }

  if (save) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = file, plot = combined, dpi = dpi,
                    width = width, height = height, units = "in", limitsize = FALSE)
  }

  return(combined)
}

# ------------------------------------------------------------
# plotTraitDendrogramFromPC
# ------------------------------------------------------------
# Builds a trait dendrogram from the PC-trait effect-size matrix.
# Traits are clustered using their association profiles across PCs.
# ------------------------------------------------------------
plotTraitDendrogramFromPC <- function(pc_trait_stats,
                                      output_file,
                                      cor_column = c("bicor", "cor"),
                                      p_column = "p",
                                      width = 10,
                                      height = 8,
                                      verbose = TRUE) {
  cor_column <- match.arg(cor_column)

  df <- .standardize_pc_trait_table(
    pc_trait_stats = pc_trait_stats,
    cor_column = cor_column,
    p_column = p_column
  )

  effect_wide <- reshape2::dcast(df, trait ~ PC, value.var = "effect")
  rownames(effect_wide) <- effect_wide$trait
  effect_wide$trait <- NULL
  effect_mat <- as.matrix(effect_wide)

  if (nrow(effect_mat) < 2) {
    stop("Need at least 2 traits to build a dendrogram.")
  }

  hc <- stats::hclust(stats::dist(effect_mat))

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  grDevices::pdf(output_file, width = width, height = height)
  plot(
    hc,
    main = paste0("Trait Dendrogram from PC Associations (", cor_column, ")"),
    xlab = "",
    sub = ""
  )
  grDevices::dev.off()

  if (verbose) {
    message("[plotTraitDendrogramFromPC] Saved: ", output_file)
  }

  return(invisible(hc))
}

# ============================================================
# Shared helper functions for ME-trait analysis / plotting
# ============================================================

# ------------------------------------------------------------
# parse_bool
#   Coerce a CLI string (TRUE/T/1/yes/y, FALSE/F/0/no/n; any case)
#   to a logical; stops with arg_name in the message otherwise.
# ------------------------------------------------------------
parse_bool <- function(x, arg_name = "argument") {
  if (is.logical(x)) return(x)
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE. Got: ", x)
}

# ------------------------------------------------------------
# write_log_lines
#   Write a character vector to `file`, one element per line
#   (overwrites).
# ------------------------------------------------------------
write_log_lines <- function(lines, file) {
  writeLines(as.character(lines), con = file)
}

# ------------------------------------------------------------
# write_vector_file
#   Write the unique values of x to `file`, one per line (empty file
#   if x is empty). Used for module/trait lists.
# ------------------------------------------------------------
write_vector_file <- function(x, file) {
  x <- unique(as.character(x))
  if (length(x) == 0) x <- character(0)
  writeLines(x, con = file)
}

# ------------------------------------------------------------
# validate_modules_object
#   Check a comethyl module object: list with numeric $MEs
#   (samples x modules), >= 2 samples, >= 1 module, unique sample IDs in
#   rownames. Returns the object with $MEs coerced to a matrix.
# ------------------------------------------------------------
validate_modules_object <- function(x, label = "modules object") {
  if (!is.list(x)) stop(label, " must be a list-like module object.")
  if (is.null(x$MEs)) stop(label, " is missing $MEs.")
  if (!(is.matrix(x$MEs) || is.data.frame(x$MEs))) stop(label, "$MEs must be matrix-like.")

  MEs <- as.matrix(x$MEs)
  if (!is.numeric(MEs)) stop(label, "$MEs must be numeric.")
  if (nrow(MEs) < 2) stop(label, "$MEs must have at least 2 samples.")
  if (ncol(MEs) < 1) stop(label, "$MEs must have at least 1 module.")
  if (is.null(rownames(MEs))) stop(label, "$MEs must have sample IDs in rownames.")
  if (anyDuplicated(rownames(MEs))) {
    dup_ids <- unique(rownames(MEs)[duplicated(rownames(MEs))])
    stop(label, "$MEs has duplicated sample IDs. Example: ",
         paste(head(dup_ids, 10), collapse = ", "))
  }

  x$MEs <- MEs
  x
}
# ------------------------------------------------------------
# plotMEtraitCor
# ------------------------------------------------------------
# Builds a ME-trait heatmap and can adjust colors.
# Defaults are the publication settings from 7b5ab0b.
# ------------------------------------------------------------
plotMEtraitCor <- function(MEtraitCor,
                           moduleOrder = 1:length(unique(MEtraitCor$module)),
                           traitOrder  = 1:length(unique(MEtraitCor$trait)),
                           topOnly = FALSE, nTop = 15, p = 0.05,
                           label.type = c("star", "p"), label.size = 8,
                           label.nudge_y = -0.38,
                           colors = blueWhiteRed(100, gamma = 0.9), limit = NULL,
                           axis.text.size = 18, legend.position = c(1.08, 0.915),
                           legend.text.size = 18, legend.title.size = 18,
                           colColorMargins = c(-0.7,4.21,1.2,11.07),
                           save = TRUE,
                           file = "ME_Trait_Correlation_Heatmap.pdf",
                           width = 11, height = 9.5, verbose = TRUE,
                           ## Controls for the labeled color bar
                           showColorBar = TRUE,
                           showColorBarLabels = TRUE,             # add names on boxes
                           colorBarLabelPos = c("inside","below"),
                           colorBarLabelSize = 4,                 # font size of the names on the color boxes
                           colorBarLabelAngle = 90,               # 90 = vertical
                           colorBarRelHeight = 0.10,              # bar height rel to heatmap
                           syncWidths = TRUE,                     # force exact alignment
                           autoContrastLabels = TRUE              # black/white text auto
) {
  label.type <- match.arg(label.type)
  colorBarLabelPos <- match.arg(colorBarLabelPos)

  if (isTRUE(verbose)) {
    message("[plotMEtraitCor] Plotting ME trait correlation heatmap")
  }

  ## Relevel factors by requested order
  MEtraitCor$module <- factor(MEtraitCor$module,
                              levels = levels(MEtraitCor$module)[moduleOrder])
  MEtraitCor$trait  <- factor(MEtraitCor$trait,
                              levels = levels(MEtraitCor$trait)[rev(traitOrder)])

  ## Significant flag (no pipes for portability)
  MEtraitCor$significant <- factor((MEtraitCor$p < p) & !is.na(MEtraitCor$p),
                                   levels = c(TRUE, FALSE))

  ## Optional top-only filtering
  if (isTRUE(topOnly)) {
    if (nrow(MEtraitCor) > nTop) {
      top_mask <- MEtraitCor$p %in% sort(MEtraitCor$p)[1:nTop] & !is.na(MEtraitCor$p)
    } else {
      top_mask <- rep(TRUE, nrow(MEtraitCor))
    }
    topModules <- unique(as.character(MEtraitCor$module[top_mask]))
    topTraits  <- unique(as.character(MEtraitCor$trait [top_mask]))
    MEtraitCor <- subset(MEtraitCor, module %in% topModules & trait %in% topTraits)
    MEtraitCor$module <- factor(MEtraitCor$module,
                                levels = levels(MEtraitCor$module)[levels(MEtraitCor$module) %in% topModules])
  }

  ## Choose correlation column
  if ("bicor" %in% colnames(MEtraitCor)) {
    corType <- "bicor"
  } else if ("cor" %in% colnames(MEtraitCor)) {
    corType <- "cor"
  } else {
    stop("[plotMEtraitCor] corType unknown, must be either 'bicor' or 'cor'")
  }

  corData <- MEtraitCor[[corType]]
  if (is.null(limit)) limit <- max(abs(corData), na.rm = TRUE)
  titleName <- paste0(toupper(substring(corType, 1, 1)), substring(corType, 2))

  ## Are module names valid R colors? (use fully qualified to avoid shadowing)
  uniqMods   <- levels(MEtraitCor$module)
  colModules <- all(uniqMods %in% grDevices::colors())

  ## Extra bottom margin for heatmap when we add a color bar
  hmMarginB <- ifelse(colModules && showColorBar, yes = 2, no = 0)

  ## Legend positioning (new ggplot2 style, fallback for older versions)
  legendTheme <- if (is.numeric(legend.position)) {
    tryCatch(
      ggplot2::theme(legend.position.inside = legend.position),
      error = function(e) ggplot2::theme(legend.position = legend.position)
    )
  } else {
    ggplot2::theme(legend.position = legend.position)
  }

  ## Base heatmap
  heatmap <- ggplot2::ggplot(MEtraitCor) +
    ggplot2::geom_tile(ggplot2::aes(x = module, y = trait, color = corData, fill = corData)) +
    ggplot2::scale_fill_gradientn(titleName, colors = colors,
                                  limits = c(-limit, limit),
                                  aesthetics = c("color", "fill")) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_alpha_manual(breaks = c(TRUE, FALSE),
                                values = c(`TRUE` = 1, `FALSE` = 0),
                                guide = "none") +
    ggplot2::theme_bw(base_size = 24) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = axis.text.size, color = "black",
                                          angle = 90, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = axis.text.size, color = "black"),
      axis.ticks  = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.title  = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = legend.text.size),
      legend.title = ggplot2::element_text(size = legend.title.size),
      panel.background = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(color = "black", linewidth = 1.25),
      panel.grid       = ggplot2::element_blank(),
      plot.background  = ggplot2::element_blank(),
      plot.margin      = grid::unit(c(1, 6, hmMarginB, 1), "lines")
    ) + legendTheme

  if (label.type == "p") {
    heatmap <- heatmap +
      ggplot2::geom_text(
        ggplot2::aes(x = module, y = trait, alpha = significant,
                     label = format(p, digits = 1, scientific = TRUE)),
        color = "black", size = label.size, nudge_y = label.nudge_y
      )
  } else {
    heatmap <- heatmap +
      ggplot2::geom_text(
        ggplot2::aes(x = module, y = trait, alpha = significant),
        label = "*", color = "black", size = label.size, nudge_y = label.nudge_y
      )
  }

  # Hide x-axis labels whenever the color bar is requested (7b5ab0b behavior;
  # originally `showColorBar && showColorBarLabels`, changed so the bar replaces
  # the axis labels even when names are not drawn on the boxes).
  if (showColorBar) {
    heatmap <- heatmap +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  ## Optional labeled color bar (boxes + names) aligned to columns
  colColors <- NULL
  if (colModules && showColorBar) {
    if (isTRUE(verbose)) {
      message("[plotMEtraitCor] Drawing module color bar",
              if (showColorBarLabels) " with labels" else "")
    }

    dfBar <- data.frame(
      module = factor(uniqMods, levels = uniqMods),
      y = 1
    )

    ## auto-contrast label color for readability
    if (showColorBarLabels) {
      if (isTRUE(autoContrastLabels)) {
        lum <- function(cols) {
          rgb <- grDevices::col2rgb(cols) / 255
          drop(0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ])
        }
        dfBar$labcol <- ifelse(lum(as.character(dfBar$module)) > 0.5, "black", "white")
      } else {
        dfBar$labcol <- "black"
      }
    }

    baseBar <- ggplot2::ggplot(dfBar, ggplot2::aes(x = module, y = y, fill = module)) +
      ggplot2::geom_tile(width = 0.98, height = 0.8) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = grid::unit(colColorMargins, "lines")
      )

    if (isTRUE(showColorBarLabels)) {
      if (colorBarLabelPos == "inside") {
        baseBar <- baseBar +
          ggplot2::geom_text(
            ggplot2::aes(label = as.character(module), colour = labcol),
            angle = colorBarLabelAngle, vjust = 0.5, size = colorBarLabelSize
          ) +
          ggplot2::scale_colour_identity()
      } else { # "below"
        baseBar <- baseBar +
          ggplot2::geom_text(
            ggplot2::aes(y = 0.05, label = as.character(module), colour = labcol),
            angle = 0, vjust = 1, size = colorBarLabelSize
          ) +
          ggplot2::scale_colour_identity() +
          ggplot2::ylim(0, 1.35)
      }
    }

    colColors <- baseBar
  }

  ## Combine and sync widths for exact alignment
  if (!is.null(colColors)) {
    if (isTRUE(syncWidths)) {
      g_heat  <- ggplot2::ggplotGrob(heatmap)
      g_strip <- ggplot2::ggplotGrob(colColors)
      g_strip$widths <- g_heat$widths
      combined <- gridExtra::arrangeGrob(
        g_heat, g_strip, ncol = 1,
        heights = grid::unit.c(grid::unit(1, "null"),
                               grid::unit(colorBarRelHeight, "null"))
      )
    } else {
      combined <- cowplot::plot_grid(heatmap, colColors, nrow = 2,
                                     rel_heights = c(1, colorBarRelHeight))
    }
  } else {
    combined <- heatmap
  }

  if (isTRUE(save)) {
    if (isTRUE(verbose)) {
      message("[plotMEtraitCor] Saving plot as ", file)
    }
    ggplot2::ggsave(filename = file, plot = combined, dpi = 600,
                    width = width, height = height, units = "in")
  }

  return(combined)
}

# ------------------------------------------------------------
# save_me_trait_method_outputs
#   Write all outputs for one ME-trait correlation method
#   (Pearson or Bicor) to out_dir: full / significant (p < p_thresh) /
#   outcome-only TSVs, a multi-sheet xlsx, and FULL + TOP heatmaps via
#   plotMEtraitCor() (columns ordered by moduleDendro$order).
#   Returns the four tables invisibly.
# ------------------------------------------------------------
save_me_trait_method_outputs <- function(MEtraitCor,
                                         method_name,
                                         out_dir,
                                         moduleDendro,
                                         p_thresh,
                                         top_n,
                                         outcome_traits_found = character(0)) {

  method_name <- match.arg(method_name, c("Pearson", "Bicor"))

  full_tsv <- file.path(out_dir, paste0("ME_Trait_Correlation_Stats_", method_name, ".tsv"))
  sig_tsv <- file.path(out_dir, paste0("ME_Trait_Correlation_Stats_", method_name, "_significant.tsv"))
  outcome_tsv <- file.path(out_dir, paste0("ME_Trait_Correlation_Stats_", method_name, "_outcome_only.tsv"))
  outcome_sig_tsv <- file.path(out_dir, paste0("ME_Trait_Correlation_Stats_", method_name, "_outcome_only_significant.tsv"))
  xlsx_file <- file.path(out_dir, paste0("ME_Trait_Correlations_", method_name, ".xlsx"))
  full_pdf <- file.path(out_dir, paste0("ME_Trait_Correlation_Heatmap_", method_name, "_FULL.pdf"))
  top_pdf <- file.path(out_dir, paste0("ME_Trait_Correlation_Heatmap_", method_name, "_TOP.pdf"))

  readr::write_tsv(MEtraitCor, full_tsv)

  sig_df <- MEtraitCor[!is.na(MEtraitCor$p) & MEtraitCor$p < p_thresh, , drop = FALSE]
  readr::write_tsv(sig_df, sig_tsv)

  outcome_df <- MEtraitCor[MEtraitCor$trait %in% outcome_traits_found, , drop = FALSE]
  outcome_sig_df <- outcome_df[!is.na(outcome_df$p) & outcome_df$p < p_thresh, , drop = FALSE]

  if (length(outcome_traits_found) > 0) {
    readr::write_tsv(outcome_df, outcome_tsv)
    readr::write_tsv(outcome_sig_df, outcome_sig_tsv)
  }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "All_Correlations")
  openxlsx::writeData(wb, "All_Correlations", MEtraitCor)

  openxlsx::addWorksheet(wb, "Significant")
  openxlsx::writeData(wb, "Significant", sig_df)

  if (length(outcome_traits_found) > 0) {
    openxlsx::addWorksheet(wb, "Outcome_Only")
    openxlsx::writeData(wb, "Outcome_Only", outcome_df)

    openxlsx::addWorksheet(wb, "Outcome_Only_Significant")
    openxlsx::writeData(wb, "Outcome_Only_Significant", outcome_sig_df)
  }

  openxlsx::saveWorkbook(wb, file = xlsx_file, overwrite = TRUE)

  plotMEtraitCor(
    MEtraitCor,
    moduleOrder = moduleDendro$order,
    p = p_thresh,
    topOnly = FALSE,
    file = full_pdf,
    width = 11,
    height = 9.5,
    colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
  )

  plotMEtraitCor(
    MEtraitCor,
    moduleOrder = moduleDendro$order,
    topOnly = TRUE,
    nTop = top_n,
    p = p_thresh,
    label.type = "p",
    label.size = 4,
    label.nudge_y = 0,
    legend.position = c(1.11, 0.795),
    colColorMargins = c(-1, 4.75, 0.5, 10.1),
    file = top_pdf,
    width = 8.5,
    height = 4.25
  )

  invisible(list(
    full = MEtraitCor,
    significant = sig_df,
    outcome = outcome_df,
    outcome_significant = outcome_sig_df
  ))
}

# ============================================================
# Shared helpers for ME-trait pair plots
# ============================================================

# ------------------------------------------------------------
# auto_ylim
#   y-limits for an eigengene plot: finite range of ME_vec padded by
#   8% (falls back to +/-0.1 for constant input, c(-1, 1) if empty).
# ------------------------------------------------------------
auto_ylim <- function(ME_vec) {
  y <- ME_vec[is.finite(ME_vec)]
  if (length(y) == 0) return(c(-1, 1))
  r <- range(y)
  pad <- 0.08 * diff(r)
  if (!is.finite(pad) || pad == 0) pad <- 0.1
  c(r[1] - pad, r[2] + pad)
}

# ------------------------------------------------------------
# is_binary_like
#   TRUE if x has exactly two distinct non-NA values.
# ------------------------------------------------------------
is_binary_like <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) == 2
}

# ------------------------------------------------------------
# relabel_binary_to_yesno
#   Recode a two-level variable to factor No/Yes. 0/1 and 1/2 map
#   low -> No; logical FALSE -> No; otherwise the sorted levels map to
#   No/Yes. Non-binary input is returned unchanged.
# ------------------------------------------------------------
relabel_binary_to_yesno <- function(x) {
  if (is.factor(x)) {
    lev <- levels(x)
    if (length(lev) != 2) return(x)
    if (all(lev %in% c("0", "1"))) return(factor(x, levels = c("0", "1"), labels = c("No", "Yes")))
    if (all(lev %in% c("1", "2"))) return(factor(x, levels = c("1", "2"), labels = c("No", "Yes")))
    if (all(tolower(lev) %in% c("false", "true"))) {
      return(factor(x, levels = lev, labels = c("No", "Yes")))
    }
    return(factor(x, levels = lev, labels = c("No", "Yes")))
  }

  ux <- unique(x[!is.na(x)])
  if (length(ux) != 2) return(x)

  if (is.numeric(x) || is.integer(x)) {
    if (all(sort(ux) == c(0, 1))) return(factor(x, levels = c(0, 1), labels = c("No", "Yes")))
    if (all(sort(ux) == c(1, 2))) return(factor(x, levels = c(1, 2), labels = c("No", "Yes")))
    uxs <- sort(ux)
    return(factor(x, levels = uxs, labels = c("No", "Yes")))
  }

  if (is.logical(x)) {
    return(factor(x, levels = c(FALSE, TRUE), labels = c("No", "Yes")))
  }

  uxs <- sort(as.character(ux))
  factor(as.character(x), levels = uxs, labels = c("No", "Yes"))
}

# ------------------------------------------------------------
# make_me_trait_subtitle
#   "r = <cor>, p = <p>" subtitle for one module x trait pair,
#   looked up in MEtraitCor (bicor preferred over cor). NULL if the
#   pair is not found exactly once or values are not finite.
# ------------------------------------------------------------
make_me_trait_subtitle <- function(MEtraitCor, module, trait) {
  if (is.null(MEtraitCor) || nrow(MEtraitCor) == 0) return(NULL)

  hit <- MEtraitCor[MEtraitCor$module == module & MEtraitCor$trait == trait, , drop = FALSE]
  if (nrow(hit) != 1) return(NULL)

  cor_col <- if ("bicor" %in% names(hit)) {
    "bicor"
  } else if ("cor" %in% names(hit)) {
    "cor"
  } else {
    NA_character_
  }

  if (is.na(cor_col)) return(NULL)

  r <- suppressWarnings(as.numeric(hit[[cor_col]][1]))
  p <- suppressWarnings(as.numeric(hit$p[1]))

  if (!is.finite(r) || !is.finite(p)) return(NULL)

  p_txt <- ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2)))
  paste0("r = ", signif(r, 3), ", ", p_txt)
}

# ------------------------------------------------------------
# add_plot_subtitle
#   Add a centred italic subtitle to a ggplot; returns p unchanged
#   when subtitle_text or p is NULL.
# ------------------------------------------------------------
add_plot_subtitle <- function(p, subtitle_text) {
  if (is.null(subtitle_text) || is.null(p)) return(p)

  p +
    ggplot2::labs(subtitle = subtitle_text) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        face = "italic",
        margin = ggplot2::margin(b = 10)
      )
    )
}

# ------------------------------------------------------------
# plotMEtraitViolin
#   Violin + boxplot + jitter of module eigengene values by a
#   categorical trait (No/Yes coloured no_color/yes_color). Drops NA
#   pairs; returns NULL if nothing is left. Does not save.
# ------------------------------------------------------------
plotMEtraitViolin <- function(ME_vec,
                              trait_fac,
                              module_name,
                              trait_name,
                              ylim,
                              no_color = "blue",
                              yes_color = "red",
                              base_size = 16) {

  df <- data.frame(
    ME = as.numeric(ME_vec),
    trait = trait_fac
  ) %>%
    dplyr::filter(!is.na(trait), !is.na(ME))

  if (nrow(df) == 0) return(NULL)

  if (!is.factor(df$trait)) {
    df$trait <- as.factor(df$trait)
  }

  if (all(levels(df$trait) %in% c("No", "Yes"))) {
    df$trait <- factor(df$trait, levels = c("No", "Yes"))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = trait, y = ME, fill = trait, color = trait)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.55, linewidth = 0.8) +
    ggplot2::geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.6
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      height = 0,
      color = "black",
      size = 1.3,
      alpha = 0.85
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = trait_name,
      y = paste0(module_name, " Module Eigengene")
    ) +
    ggplot2::coord_cartesian(ylim = ylim)

  levs <- levels(df$trait)
  if (length(levs) == 2 && all(c("No", "Yes") %in% levs)) {
    p <- p +
      ggplot2::scale_fill_manual(values = c("No" = no_color, "Yes" = yes_color), drop = FALSE) +
      ggplot2::scale_color_manual(values = c("No" = no_color, "Yes" = yes_color), drop = FALSE)
  }

  p
}

# ------------------------------------------------------------
# collect_set_files
#   Gather trait-set .txt files from up to five explicit paths and/or
#   every .txt in set_dir. Stops on missing files or an empty set.
#   Returns unique normalized paths.
# ------------------------------------------------------------
collect_set_files <- function(set_dir = NULL,
                              set_file = NULL,
                              set_file2 = NULL,
                              set_file3 = NULL,
                              set_file4 = NULL,
                              set_file5 = NULL) {

  direct_files <- c(set_file, set_file2, set_file3, set_file4, set_file5)
  direct_files <- direct_files[!is.na(direct_files) & nzchar(direct_files)]

  if (length(direct_files) > 0) {
    missing_direct <- direct_files[!file.exists(direct_files)]
    if (length(missing_direct) > 0) {
      stop("The following trait set files do not exist:\n  ",
           paste(missing_direct, collapse = "\n  "))
    }
  }

  dir_files <- character(0)
  if (!is.null(set_dir) && nzchar(set_dir)) {
    if (!dir.exists(set_dir)) {
      stop("set_dir does not exist: ", set_dir)
    }

    dir_files <- list.files(
      set_dir,
      pattern = "\\.txt$",
      full.names = TRUE
    )

    if (length(dir_files) == 0) {
      stop("No .txt files found in set_dir: ", set_dir)
    }
  }

  all_files <- c(direct_files, dir_files)

  if (length(all_files) == 0) {
    stop("No trait set files provided. Use --set_dir or at least one --set_file.")
  }

  unique(normalizePath(all_files, mustWork = TRUE))
}

# ------------------------------------------------------------
# read_trait_set_file
#   Read one trait-set file: trimmed, non-blank, unique lines.
# ------------------------------------------------------------
read_trait_set_file <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

# ------------------------------------------------------------
# get_set_name
#   Set name = file name without directory or extension.
# ------------------------------------------------------------
get_set_name <- function(file) {
  tools::file_path_sans_ext(basename(file))
}

# ================================================================
# ANNOTATION HELPERS -- CLI / LOGGING / VALIDATION + LEGACY ANNOTATION PATH
#
# Gene assignment in the CURRENT pipeline is GREAT (see the "GREAT REGION
# -> GENE" section further down): 12a --gene_assignment great (default)
# gives up to two genes per region with comethyl::annotateModule()'s
# basalPlusExt rule, run offline. --gene_assignment nearest instead uses
# offline_nearest_gene() below (one gene per region, nearest gene body).
#
# annotate_offline_only() / annotate_regions_safe() / ChIPseeker
# genomic_location / DMRichR-style CpG flags are the PRE-annotatr 12a
# path. The current 12a does not call them; they are kept so earlier
# outputs (and any script still using them, e.g. an older
# 12c_annotate_region_matrix.R) remain reproducible.
#
# All Bioconductor calls use pkg::fn() namespacing, so those packages
# only need to be installed, not attached, in the calling script.
# ================================================================

# ================================================================
# REPRODUCIBILITY HELPERS
#
# write_session_info(): captures exact R + package versions and the
#   git commit hash of the pipeline code, so a given run_log.txt can
#   be tied to the exact code + environment that produced it.
# file_checksum(): md5 of an input file, so you can prove downstream
#   which exact copy of e.g. Filtered_Regions.txt or cpgIslandExt.txt
#   was used, independent of the path (which can drift/get renamed).
# ================================================================

write_session_info <- function(file, extra_files = NULL) {
  git_hash <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  if (length(git_hash) == 0 || !nzchar(git_hash)) git_hash <- "(not a git repo or git unavailable)"

  git_dirty <- tryCatch({
    status <- system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)
    if (length(status) == 0) "clean" else "DIRTY (uncommitted changes present)"
  }, error = function(e) "(unknown)")

  lines <- c(
    paste0("timestamp\t", timestamp_now()),
    paste0("git_commit\t", git_hash[1]),
    paste0("git_status\t", git_dirty),
    ""
  )

  if (!is.null(extra_files)) {
    lines <- c(lines, "# Input file checksums (md5)", "")
    for (nm in names(extra_files)) {
      f <- extra_files[[nm]]
      if (!is.null(f) && file.exists(f)) {
        lines <- c(lines, paste0(nm, "\t", f, "\t", tools::md5sum(f)))
      } else {
        lines <- c(lines, paste0(nm, "\t", ifelse(is.null(f), "(not provided)", f), "\t(file not found)"))
      }
    }
    lines <- c(lines, "")
  }

  lines <- c(lines, "# sessionInfo()", "", utils::capture.output(print(sessionInfo())))

  write_lines_safe(lines, file)
}

# ------------------------------------------------------------
# file_checksum
#   md5 of a file, or NA if path is NULL / missing.
# ------------------------------------------------------------
file_checksum <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

# ------------------------------------------------------------
# CLI argument helpers
# ------------------------------------------------------------
# ------------------------------------------------------------
# get_arg
#   Value following `flag` in commandArgs(trailingOnly = TRUE), else
#   default. (Scripts using optparse don't need this.)
# ------------------------------------------------------------
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

# ------------------------------------------------------------
# trim_or_null
#   Trim whitespace; NULL, NA and "" become NULL.
# ------------------------------------------------------------
trim_or_null <- function(x) {
  if (is.null(x) || is.na(x)) return(NULL)
  x <- trimws(x)
  if (!nzchar(x)) return(NULL)
  x
}

# ------------------------------------------------------------
# Filesystem / logging helpers
# ------------------------------------------------------------
# ------------------------------------------------------------
# safe_dir_create
#   Create a directory (recursively) if it doesn't exist.
# ------------------------------------------------------------
safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# timestamp_now
#   Current time as "YYYY-mm-dd HH:MM:SS".
# ------------------------------------------------------------
timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

# ------------------------------------------------------------
# append_log
#   Paste ... into one timestamped line, print it, and append it to
#   logfile.
# ------------------------------------------------------------
append_log <- function(logfile, ...) {
  txt <- paste0("[", timestamp_now(), "] ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = logfile, append = TRUE)
}

# ------------------------------------------------------------
# write_lines_safe
#   writeLines() with useBytes = TRUE (overwrites `file`).
# ------------------------------------------------------------
write_lines_safe <- function(x, file) {
  writeLines(as.character(x), con = file, useBytes = TRUE)
}

# ------------------------------------------------------------
# Validation helpers
# ------------------------------------------------------------
# ------------------------------------------------------------
# stop_if_missing
#   Stop with "Missing required argument: <label>" if x is NULL or "".
# ------------------------------------------------------------
stop_if_missing <- function(x, label) {
  if (is.null(x) || !nzchar(x)) stop("Missing required argument: ", label, call. = FALSE)
}

# ------------------------------------------------------------
# validate_file_exists
#   Stop with "<label> not found: <path>" if path doesn't exist.
# ------------------------------------------------------------
validate_file_exists <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

# req_cols lets callers relax the requirement (e.g. drop "module"
# for matrix-level annotation where there is no module assignment yet)
validate_regions_df <- function(regions, source_label = "regions",
                                 req_cols = c("RegionID", "chr", "start", "end", "module")) {
  missing_cols <- setdiff(req_cols, colnames(regions))
  if (length(missing_cols) > 0) {
    stop(
      source_label, " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

# ------------------------------------------------------------
# Extract a regions data.frame (RegionID, chr, start, end, module)
# from a comethyl module RDS object -- shared by 12a (module-level
# annotation) and 12c (merging module assignment into per-variant
# region x sample tables).
# ------------------------------------------------------------
extract_regions_from_module_object <- function(obj, label = "module object") {
  if (is.list(obj) && "regions" %in% names(obj)) {
    regions <- obj$regions
  } else if (is.data.frame(obj)) {
    regions <- obj
  } else {
    stop(
      "Could not extract regions from ", label,
      ". Expected either a list with $regions or a data.frame.",
      call. = FALSE
    )
  }

  validate_regions_df(regions, label)
  regions
}

# ------------------------------------------------------------
# Genome -> annotation package resolution
#
# Mirrors DMRichR::annotationDatabases()'s genome -> TxDb/org.db
# mapping, so --genome actually determines which packages get used
# instead of them being hardcoded to hg38/human. Two deliberate
# differences from DMRichR's version:
#   - No BSgenome/sequence loading -- nothing here needs sequence,
#     only transcript models (TxDb) and gene ID mappings (org.db).
#   - No live BiocManager::install() of missing packages. A SLURM
#     job silently installing packages mid-run is exactly the kind
#     of non-reproducible, network-dependent behavior this pipeline
#     has been fixed to avoid (see the GREAT/rGREAT saga). If a
#     required package isn't installed, this errors with a clear
#     message telling you what to add to the pixi environment,
#     rather than trying to fetch it at runtime.
# ------------------------------------------------------------
resolve_annotation_packages <- function(genome) {
  mapping <- list(
    hg38     = list(txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",      annoDb = "org.Hs.eg.db"),
    hg19     = list(txdb = "TxDb.Hsapiens.UCSC.hg19.knownGene",      annoDb = "org.Hs.eg.db"),
    mm10     = list(txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene",     annoDb = "org.Mm.eg.db"),
    mm9      = list(txdb = "TxDb.Mmusculus.UCSC.mm9.knownGene",      annoDb = "org.Mm.eg.db"),
    rheMac10 = list(txdb = "TxDb.Mmulatta.UCSC.rheMac10.refGene",    annoDb = "org.Mmu.eg.db"),
    rheMac8  = list(txdb = "TxDb.Mmulatta.UCSC.rheMac8.refGene",     annoDb = "org.Mmu.eg.db"),
    rn6      = list(txdb = "TxDb.Rnorvegicus.UCSC.rn6.refGene",      annoDb = "org.Rn.eg.db"),
    danRer11 = list(txdb = "TxDb.Drerio.UCSC.danRer11.refGene",      annoDb = "org.Dr.eg.db"),
    galGal6  = list(txdb = "TxDb.Ggallus.UCSC.galGal6.refGene",      annoDb = "org.Gg.eg.db"),
    bosTau9  = list(txdb = "TxDb.Btaurus.UCSC.bosTau9.refGene",      annoDb = "org.Bt.eg.db"),
    panTro6  = list(txdb = "TxDb.Ptroglodytes.UCSC.panTro6.refGene", annoDb = "org.Pt.eg.db"),
    dm6      = list(txdb = "TxDb.Dmelanogaster.UCSC.dm6.ensGene",    annoDb = "org.Dm.eg.db"),
    susScr11 = list(txdb = "TxDb.Sscrofa.UCSC.susScr11.refGene",     annoDb = "org.Ss.eg.db"),
    canFam3  = list(txdb = "TxDb.Cfamiliaris.UCSC.canFam3.refGene",  annoDb = "org.Cf.eg.db"),
    TAIR10   = list(txdb = "TxDb.Athaliana.BioMart.plantsmart28",    annoDb = "org.At.tair.db"),
    TAIR9    = list(txdb = "TxDb.Athaliana.BioMart.plantsmart28",    annoDb = "org.At.tair.db")
  )

  if (!genome %in% names(mapping)) {
    stop(
      "genome '", genome, "' is not supported. Choose one of: ",
      paste(names(mapping), collapse = ", "), " [case sensitive].",
      call. = FALSE
    )
  }

  mapping[[genome]]
}

# ------------------------------------------------------------
# offline_nearest_gene
# ------------------------------------------------------------
# Single nearest gene per region (12a --gene_assignment nearest; also the
# legacy offline fallback). NOT comethyl's GREAT rule.
#   - Genes = whole gene-body spans (union of transcripts), not TSSs.
#   - distanceToNearest(ignore.strand = TRUE): one gene per region; a
#     region inside a gene always gets that gene; ties/overlapping genes
#     are broken arbitrarily; no maximum distance.
#   - hg38 uses EnsDb.Hsapiens.v86 if installed, else TxDb knownGene.
#     Symbols come from org.Hs.eg.db either way (NA when an Ensembl gene
#     has no Entrez mapping). CAUTION: EnsDb uses "1"-style seqnames;
#     with "chr1"-style regions that branch finds no chromosomes in common.
# Returns data.frame(chr, start, end, gene_symbol, gene_entrezID).
# ------------------------------------------------------------
offline_nearest_gene <- function(gr, genome = "hg38", verbose = TRUE) {
  pkgs <- resolve_annotation_packages(genome)
  txdb_pkg  <- pkgs$txdb
  annoDb_pkg <- pkgs$annoDb

  if (!requireNamespace(annoDb_pkg, quietly = TRUE)) {
    stop(annoDb_pkg, " is required for offline annotation of genome '", genome,
         "'. Install it in the pixi environment.", call. = FALSE)
  }

  # hg38 keeps its existing EnsDb.Hsapiens.v86-preferred path unchanged
  # (matches prior validated behavior exactly). Every other genome goes
  # straight through the genome-resolved TxDb + org.db path.
  use_ensdb <- genome == "hg38" && requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)

  if (use_ensdb) {
    if (verbose) message("[offline] Using EnsDb.Hsapiens.v86")
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    seqs <- unique(as.character(GenomeInfoDb::seqnames(gr)))
    seqs <- seqs[seqs %in% GenomeInfoDb::seqlevels(edb)]

    genes <- ensembldb::genes(edb, filter = AnnotationFilter::SeqNameFilter(seqs))
    ggr <- GenomicRanges::GRanges(genes)

    ens_ids <- genes$gene_id
    annoDb_obj <- get(annoDb_pkg, envir = asNamespace(annoDb_pkg))
    map <- AnnotationDbi::select(
      annoDb_obj,
      keys = ens_ids,
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "ENTREZID")
    )

    ggr$SYMBOL   <- map$SYMBOL[match(genes$gene_id, map$ENSEMBL)]
    ggr$ENTREZID <- map$ENTREZID[match(genes$gene_id, map$ENSEMBL)]

  } else {
    if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
      stop(txdb_pkg, " is required for offline annotation of genome '", genome,
           "'. Install it in the pixi environment.", call. = FALSE)
    }
    if (verbose) message("[offline] Using ", txdb_pkg)
    txdb <- get(txdb_pkg, envir = asNamespace(txdb_pkg))
    ggr  <- GenomicFeatures::genes(txdb)

    annoDb_obj <- get(annoDb_pkg, envir = asNamespace(annoDb_pkg))
    map <- AnnotationDbi::select(
      annoDb_obj,
      keys = ggr$gene_id,
      keytype = "ENTREZID",
      columns = c("SYMBOL")
    )

    ggr$ENTREZID <- ggr$gene_id
    ggr$SYMBOL   <- map$SYMBOL[match(ggr$ENTREZID, map$ENTREZID)]
  }

  hit <- GenomicRanges::distanceToNearest(gr, ggr, ignore.strand = TRUE)
  ng  <- ggr[S4Vectors::subjectHits(hit)]

  data.frame(
    chr = as.character(GenomeInfoDb::seqnames(gr))[S4Vectors::queryHits(hit)],
    start = as.integer(S4Vectors::start(gr))[S4Vectors::queryHits(hit)],
    end = as.integer(S4Vectors::end(gr))[S4Vectors::queryHits(hit)],
    gene_symbol = as.character(ng$SYMBOL),
    gene_entrezID = as.character(ng$ENTREZID),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# CpG island context (island / shore / shelf / open sea) --
# replicates DMRichR::getCpGs()'s exact algorithm (island -> stretch
# 4000bp -> setdiff for shores -> stretch again -> setdiff for
# shelves -> gaps for open sea) but reads the CpG island track from
# a LOCAL cached file instead of a live download every call, so
# results are fixed and reproducible.
#
# cpg_island_file: local UCSC cpgIslandExt table. Download once from
# ANY machine with internet (your laptop is fine -- no cluster
# internet required) and scp it over:
#
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz
#
# Then point --cpg_island_file at the local (gzipped is fine) copy.
# Uses only base GenomicRanges (resize/setdiff/gaps) -- no plyranges
# dependency, matching the "replicate locally" choice.
# ------------------------------------------------------------
build_cpg_annotations <- function(cpg_island_file, shore_bp = 2000, shelf_bp = 2000, verbose = TRUE) {
  validate_file_exists(cpg_island_file, "cpg_island_file")
  cols <- c("bin", "chrom", "chromStart", "chromEnd", "name", "length",
            "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
  df <- read.table(cpg_island_file, sep = "\t", stringsAsFactors = FALSE,
                    col.names = cols, quote = "")
  if (verbose) message("[cpg] Loaded ", nrow(df), " CpG islands from: ", cpg_island_file)

  # UCSC chromStart is 0-based half-open; GRanges is 1-based inclusive
  islands <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$chromStart + 1, end = df$chromEnd)
  )
  islands <- GenomeInfoDb::keepStandardChromosomes(islands, pruning.mode = "coarse")

  if (verbose) message("[cpg] Building CpG shores...")
  # resize(..., fix="center") with width += 2*shore_bp is equivalent to
  # DMRichR's plyranges::stretch(shore_bp*2) (extends shore_bp on each side)
  shores_raw <- GenomicRanges::resize(
    islands, width = GenomicRanges::width(islands) + 2 * shore_bp, fix = "center"
  )
  shores_raw <- GenomicRanges::trim(shores_raw)
  shores <- GenomicRanges::setdiff(shores_raw, islands)

  if (verbose) message("[cpg] Building CpG shelves...")
  shelves_raw <- GenomicRanges::resize(
    shores, width = GenomicRanges::width(shores) + 2 * shelf_bp, fix = "center"
  )
  shelves_raw <- GenomicRanges::trim(shelves_raw)
  shelves <- GenomicRanges::setdiff(shelves_raw, islands)
  shelves <- GenomicRanges::setdiff(shelves, shores)

  if (verbose) message("[cpg] Building open sea (inter-CGI)...")
  inter_cgi <- GenomicRanges::gaps(GenomicRanges::sort(c(islands, shores, shelves)))
  # gaps() also returns the space before the first range / after the last on
  # each strand-agnostic "*" track; keep only ranges actually on real chroms
  inter_cgi <- inter_cgi[GenomeInfoDb::seqnames(inter_cgi) %in% GenomeInfoDb::seqnames(islands)]

  islands$type <- "islands"
  shores$type  <- "shores"
  shelves$type <- "shelves"
  inter_cgi$type <- "inter"

  cpgs <- GenomicRanges::sort(c(islands, shores, shelves, inter_cgi))

  if (verbose) {
    message("[cpg] islands=", length(islands), " shores=", length(shores),
            " shelves=", length(shelves), " inter=", length(inter_cgi))
  }

  cpgs
}

# Four independent Yes/No overlap flags, matching DMRichR's
# CpG.Island / CpG.Shore / CpG.Shelf / Open.Sea columns exactly
# (a region can flag "Yes" in more than one if it spans categories).
classify_cpg_flags <- function(gr, cpgs, verbose = TRUE) {
  islands <- cpgs[cpgs$type == "islands"]
  shores  <- cpgs[cpgs$type == "shores"]
  shelves <- cpgs[cpgs$type == "shelves"]
  inter   <- cpgs[cpgs$type == "inter"]

  flags <- data.frame(
    CpG.Island = ifelse(GenomicRanges::countOverlaps(gr, islands) > 0, "Yes", "No"),
    CpG.Shore  = ifelse(GenomicRanges::countOverlaps(gr, shores)  > 0, "Yes", "No"),
    CpG.Shelf  = ifelse(GenomicRanges::countOverlaps(gr, shelves) > 0, "Yes", "No"),
    Open.Sea   = ifelse(GenomicRanges::countOverlaps(gr, inter)   > 0, "Yes", "No"),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message("[cpg] Island=", sum(flags$CpG.Island == "Yes"),
            " Shore=", sum(flags$CpG.Shore == "Yes"),
            " Shelf=", sum(flags$CpG.Shelf == "Yes"),
            " OpenSea=", sum(flags$Open.Sea == "Yes"))
  }

  flags
}

# Convenience single-label column derived from the 4 flags above,
# with an explicit, documented precedence rule (Island > Shore >
# Shelf > OpenSea) for anywhere you want one category per region
# (e.g. a single heatmap axis). The 4 flags above remain the
# authoritative, DMRichR-matching columns.
derive_cpg_context_label <- function(flags_df) {
  with(flags_df, ifelse(CpG.Island == "Yes", "Island",
                  ifelse(CpG.Shore  == "Yes", "Shore",
                  ifelse(CpG.Shelf  == "Yes", "Shelf", "OpenSea"))))
}

# ------------------------------------------------------------
# Genomic location (promoter / 5' UTR / exon / intron / 3' UTR /
# downstream / intergenic) -- fully local via TxDb, no internet
# needed at call time (only the one-time package install does).
# Matches DMRichR::annotateRegions()'s ChIPseeker call (overlap =
# "all") and its category-string simplification exactly, so
# genomic_location strings are byte-identical to DMRichR's own
# "annotation" column (e.g. "Promoter", "5' UTR", "Distal Intergenic").
# ------------------------------------------------------------
annotate_genomic_location <- function(gr, genome = "hg38", txdb_pkg = NULL,
                                      tss_region = c(-2000, 200), verbose = TRUE) {
  if (is.null(txdb_pkg)) txdb_pkg <- resolve_annotation_packages(genome)$txdb

  n_input <- length(gr)

  empty_result <- function() {
    data.frame(
      chr = as.character(GenomeInfoDb::seqnames(gr)),
      start = as.integer(GenomicRanges::start(gr)),
      end = as.integer(GenomicRanges::end(gr)),
      genomic_location = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
    message("[genomic_location] SKIPPED -- ChIPseeker not installed/loadable in this environment. ",
            "Install via: pixi run Rscript -e \"BiocManager::install('ChIPseeker')\"")
    return(empty_result())
  }
  if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
    message("[genomic_location] SKIPPED -- ", txdb_pkg, " not installed/loadable in this environment.")
    return(empty_result())
  }

  txdb <- get(txdb_pkg, envir = asNamespace(txdb_pkg))

  if (verbose) message("[genomic_location] Running ChIPseeker::annotatePeak(overlap = \"all\") with ", txdb_pkg)
  peakAnno <- ChIPseeker::annotatePeak(
    gr, TxDb = txdb, tssRegion = tss_region, overlap = "all", verbose = FALSE
  )
  anno_df <- as.data.frame(peakAnno)

  # DMRichR-matching simplification: strip everything after " (" --
  # same regex as DMRichR::annotateRegions()'s
  # dplyr::mutate(annotation = gsub(" \\(.*","", annotation))
  anno_df$genomic_location <- gsub(" \\(.*", "", as.character(anno_df$annotation))

  # Join back by chr/start/end -- the same proven join key
  # offline_nearest_gene() already uses successfully for gene_symbol.
  # ChIPseeker's as.data.frame() output is guaranteed to carry
  # seqnames/start/end (it IS the GRanges), unlike a custom mcol like
  # RegionID, which isn't reliably preserved under that exact name.
  out <- anno_df %>%
    dplyr::rename(chr = seqnames) %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::select(chr, start, end, genomic_location) %>%
    dplyr::distinct(chr, start, end, .keep_all = TRUE)

  if (verbose) message("[genomic_location] ChIPseeker returned ", nrow(out),
                       " unique region rows for ", n_input, " input regions")

  n_missing <- n_input - sum(
    paste(as.character(GenomeInfoDb::seqnames(gr)), GenomicRanges::start(gr), GenomicRanges::end(gr)) %in%
    paste(out$chr, out$start, out$end)
  )
  if (n_missing > 0) {
    message("[genomic_location] NOTE -- ", n_missing, " of ", n_input,
            " regions got no genomic_location (likely on a chromosome/scaffold ",
            "not present in ", txdb_pkg, "). These will have genomic_location = NA after joining.")
  }

  out
}

# ------------------------------------------------------------
# annotate_offline_only
#   LEGACY (pre-annotatr 12a). Nearest-gene annotation via
#   offline_nearest_gene(), plus optional ChIPseeker genomic_location
#   (add_genomic_location) and DMRichR-style CpG flags (cpg_island_file).
#   Joins by chr/start/end; gene_description / gene_ensemblID are NA.
#   Optionally writes file_txt. Returns regions_df with added columns.
# ------------------------------------------------------------
annotate_offline_only <- function(regions_df, genome = "hg38", file_txt = NULL, verbose = TRUE,
                                  add_genomic_location = TRUE,
                                  txdb_pkg = NULL,
                                  tss_region = c(-2000, 200),
                                  cpg_island_file = NULL,
                                  shore_bp = 2000, shelf_bp = 2000) {
  validate_regions_df(regions_df, "regions_df", req_cols = c("RegionID", "chr", "start", "end"))

  gr <- GenomicRanges::GRanges(
    seqnames = regions_df$chr,
    ranges   = IRanges::IRanges(start = regions_df$start, end = regions_df$end),
    RegionID = regions_df$RegionID
  )

  ng <- offline_nearest_gene(gr, genome = genome, verbose = verbose)

  out <- regions_df %>%
    dplyr::left_join(
      ng %>% dplyr::select(chr, start, end, gene_symbol, gene_entrezID),
      by = c("chr", "start", "end")
    ) %>%
    dplyr::mutate(
      gene_description = NA_character_,
      gene_ensemblID   = NA_character_
    )

  if (isTRUE(add_genomic_location)) {
    loc_df <- annotate_genomic_location(
      gr, genome = genome, txdb_pkg = txdb_pkg, tss_region = tss_region, verbose = verbose
    )
    out <- out %>% dplyr::left_join(loc_df, by = c("chr", "start", "end"))
    n_missing_loc <- sum(is.na(out$genomic_location))
    if (verbose) message("[genomic_location] ", nrow(out) - n_missing_loc, " of ", nrow(out),
                         " regions have a genomic_location after join")
  }

  if (!is.null(cpg_island_file)) {
    cpgs <- build_cpg_annotations(cpg_island_file, shore_bp = shore_bp, shelf_bp = shelf_bp, verbose = verbose)
    flags <- classify_cpg_flags(gr, cpgs, verbose = verbose)
    out <- cbind(out, flags)
    out$cpg_context <- derive_cpg_context_label(flags)
  } else if (verbose) {
    message("[cpg] No cpg_island_file provided -- skipping CpG.Island/Shore/Shelf/Open.Sea/cpg_context columns.")
  }

  if (!is.null(file_txt) && nzchar(file_txt)) {
    write.table(out, file = file_txt, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  out
}

# req_cols passed through so matrix-level callers (no "module" column)
# don't get rejected by validation. add_genomic_location/cpg_island_file
# etc. are forwarded to annotate_offline_only() when the offline path
# runs (either --annotation_mode offline, or as a GREAT fallback under
# "auto").
annotate_regions_safe <- function(regions_df,
                                  genome = "hg38",
                                  annotation_mode = c("auto", "great", "offline"),
                                  file_txt = NULL,
                                  verbose = TRUE,
                                  req_cols = c("RegionID", "chr", "start", "end", "module"),
                                  add_genomic_location = TRUE,
                                  txdb_pkg = NULL,
                                  tss_region = c(-2000, 200),
                                  cpg_island_file = NULL,
                                  shore_bp = 2000, shelf_bp = 2000) {
  annotation_mode <- match.arg(annotation_mode)
  validate_regions_df(regions_df, "regions_df", req_cols = req_cols)

  offline_args <- list(
    regions_df = regions_df, genome = genome, file_txt = file_txt, verbose = verbose,
    add_genomic_location = add_genomic_location, txdb_pkg = txdb_pkg,
    tss_region = tss_region, cpg_island_file = cpg_island_file,
    shore_bp = shore_bp, shelf_bp = shelf_bp
  )

  if (annotation_mode == "offline") {
    return(do.call(annotate_offline_only, offline_args))
  }

  if (annotation_mode == "great") {
    if (verbose) message("[annotate] Using comethyl::annotateModule() only")
    return(comethyl::annotateModule(regions_df, genome = genome, file = file_txt))
  }

  tryCatch(
    {
      if (verbose) message("[annotate] Trying comethyl::annotateModule()")
      comethyl::annotateModule(regions_df, genome = genome, file = file_txt)
    },
    error = function(e) {
      message("[annotate] GREAT-based annotation failed: ", conditionMessage(e))
      message("[annotate] Falling back to offline nearest-gene annotation.")
      do.call(annotate_offline_only, offline_args)
    }
  )
}

# ================================================================
# MODULE ANNOTATION ENRICHMENT (genic + CpG), per module
#
# Mirrors DMRichR::DMRichGenic() / DMRichR::DMRichCpG(), but instead
# of testing "significant DMRs" vs "background regions", each
# comethyl module is tested against all OTHER regions (2x2 Fisher's
# exact test per module x category). FDR-adjusted across all tests
# within each function call.
#
# Requires annotated_regions to already have a "module" column and
# (for genic) a "genomic_location" column, or (for CpG) the
# CpG.Island/CpG.Shore/CpG.Shelf/Open.Sea flag columns -- i.e. the
# output of annotate_offline_only()/annotate_regions_safe() with
# add_genomic_location = TRUE and cpg_island_file provided.
# ================================================================

moduleGenicEnrichment <- function(annotated_regions,
                                  module_col = "module",
                                  location_col = "genomic_location",
                                  categories = c("Promoter", "5' UTR", "Exon", "Intron",
                                                 "3' UTR", "Downstream", "Distal Intergenic"),
                                  exclude_modules = "grey",
                                  verbose = TRUE) {
  if (!location_col %in% colnames(annotated_regions)) {
    stop("annotated_regions is missing '", location_col, "' -- ",
         "re-run annotation with add_genomic_location = TRUE.", call. = FALSE)
  }

  df <- annotated_regions
  if (!is.null(exclude_modules)) df <- df[!(df[[module_col]] %in% exclude_modules), , drop = FALSE]
  modules <- sort(unique(as.character(df[[module_col]])))

  if (verbose) message("[enrich] Genic enrichment: ", length(modules), " modules x ", length(categories), " categories")

  rows <- list()
  for (mod in modules) {
    in_mod <- df[[module_col]] == mod
    for (cat in categories) {
      has_cat <- !is.na(df[[location_col]]) & grepl(cat, df[[location_col]], fixed = TRUE)
      a <- sum(in_mod & has_cat);  b <- sum(in_mod & !has_cat)
      c <- sum(!in_mod & has_cat); d <- sum(!in_mod & !has_cat)
      ft <- tryCatch(stats::fisher.test(matrix(c(a, c, b, d), nrow = 2)), error = function(e) NULL)
      if (is.null(ft)) next
      rows[[length(rows) + 1]] <- data.frame(
        module = mod, category = cat, n_in_module = a,
        OR = unname(ft$estimate), CIlower = ft$conf.int[1], CIupper = ft$conf.int[2],
        p = ft$p.value, stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    stop("moduleGenicEnrichment(): every module x category Fisher's test failed -- ",
         "check for NA/unexpected values in '", location_col, "'.", call. = FALSE)
  }

  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p, method = "fdr")
  out
}

# ------------------------------------------------------------
# moduleCpGEnrichment
#   LEGACY. Per-module 2x2 Fisher's test for each CpG flag column
#   (CpG.Island/Shore/Shelf/Open.Sea == "Yes"), module vs all other
#   non-excluded regions; BH FDR across all tests. Same output columns as
#   moduleGenicEnrichment(). Current 12a uses moduleAnnotatrEnrichment().
# ------------------------------------------------------------
moduleCpGEnrichment <- function(annotated_regions,
                                module_col = "module",
                                cpg_cols = c(CpG.Island = "Island", CpG.Shore = "Shore",
                                            CpG.Shelf = "Shelf", Open.Sea = "OpenSea"),
                                exclude_modules = "grey",
                                verbose = TRUE) {
  missing_cols <- setdiff(names(cpg_cols), colnames(annotated_regions))
  if (length(missing_cols) > 0) {
    stop("annotated_regions is missing CpG flag column(s): ", paste(missing_cols, collapse = ", "),
         " -- re-run annotation with cpg_island_file provided.", call. = FALSE)
  }

  df <- annotated_regions
  if (!is.null(exclude_modules)) df <- df[!(df[[module_col]] %in% exclude_modules), , drop = FALSE]
  modules <- sort(unique(as.character(df[[module_col]])))

  if (verbose) message("[enrich] CpG enrichment: ", length(modules), " modules x ", length(cpg_cols), " categories")

  rows <- list()
  for (mod in modules) {
    in_mod <- df[[module_col]] == mod
    for (col in names(cpg_cols)) {
      has_cat <- !is.na(df[[col]]) & df[[col]] == "Yes"
      a <- sum(in_mod & has_cat);  b <- sum(in_mod & !has_cat)
      c <- sum(!in_mod & has_cat); d <- sum(!in_mod & !has_cat)
      ft <- tryCatch(stats::fisher.test(matrix(c(a, c, b, d), nrow = 2)), error = function(e) NULL)
      if (is.null(ft)) next
      rows[[length(rows) + 1]] <- data.frame(
        module = mod, category = unname(cpg_cols[col]), n_in_module = a,
        OR = unname(ft$estimate), CIlower = ft$conf.int[1], CIupper = ft$conf.int[2],
        p = ft$p.value, stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    stop("moduleCpGEnrichment(): every module x category Fisher's test failed -- ",
         "check for NA/unexpected values in the CpG flag columns.", call. = FALSE)
  }

  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p, method = "fdr")
  out
}

# ------------------------------------------------------------
# Module x category enrichment heatmap: color = log2(OR)
# (blue = depleted, red = enriched), "*" = FDR-significant.
# Same visual language as plotMEtraitCor()/plotPCTrait() elsewhere
# in this file, so it reads consistently alongside the rest of the
# pipeline's figures.
# ------------------------------------------------------------
plotModuleAnnotationEnrichment <- function(enrichment_df,
                                           title = "Module Annotation Enrichment",
                                           p_cut = 0.05,
                                           moduleOrder = NULL,
                                           categoryOrder = NULL,
                                           colors = WGCNA::blueWhiteRed(100, gamma = 0.9),
                                           save = TRUE,
                                           file = "Module_Annotation_Enrichment.pdf",
                                           width = NULL, height = NULL,
                                           base_size = 11,
                                           axis.text.size = 8,
                                           star.size = 4,
                                           verbose = TRUE) {
  df <- enrichment_df
  df$log2OR <- log2(df$OR)
  df$log2OR[!is.finite(df$log2OR)] <- NA  # OR of 0 or Inf from small/zero counts
  df$significant <- df$fdr < p_cut & !is.na(df$fdr)

  if (is.null(moduleOrder)) moduleOrder <- sort(unique(df$module))
  if (is.null(categoryOrder)) categoryOrder <- unique(df$category)

  df$module   <- factor(df$module, levels = moduleOrder)
  df$category <- factor(df$category, levels = rev(categoryOrder))

  limit <- suppressWarnings(max(abs(df$log2OR), na.rm = TRUE))
  if (!is.finite(limit) || limit == 0) limit <- 1

  nX <- length(levels(df$module))
  nY <- length(levels(df$category))
  if (is.null(width))  width  <- max(6, min(24, nX * 0.18 + 3))
  if (is.null(height)) height <- max(3, min(10, nY * 0.35 + 2))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = module, y = category)) +
    ggplot2::geom_tile(ggplot2::aes(fill = log2OR)) +
    ggplot2::scale_fill_gradientn("log2(OR)", colors = colors, limits = c(-limit, limit), na.value = "grey85") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(significant, "*", "")), size = star.size) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = axis.text.size, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = axis.text.size),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = title)

  if (isTRUE(save)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = file, plot = p, dpi = 600,
                    width = width, height = height, units = "in", limitsize = FALSE)
    if (verbose) message("[plot] Saved: ", file)
  }

  p
}
# ================================================================
# ANNOTATR-BASED MODULE ANNOTATION (matches comethyl::annotateModule())
#
# comethyl::annotateModule() gets genic context and CpG context from
# annotatr, using exactly these annotation codes:
#   <genome>_basicgenes, <genome>_genes_intergenic,
#   <genome>_enhancers_fantom, <genome>_cpgs
# with keepStandardChromosomes() and annotate_regions(ignore.strand =
# TRUE). The functions below reproduce that, with three additions:
#
#   1. Caching. annotatr downloads the CpG island track (UCSC) and the
#      FANTOM5 enhancers (plus an hg19->hg38 liftOver chain from
#      AnnotationHub) at build time. Hive compute nodes have no
#      internet, so build once on a LOGIN node and save to an RDS
#      (--annotatr_cache); compute-node runs then just load the RDS.
#   2. Offline fallbacks. CpG features can be built from a local
#      cpgIslandExt.txt.gz using annotatr's own algorithm, and
#      enhancers can be read from a local FANTOM5 BED (+ optional
#      liftOver chain), so the cache can also be built with no internet.
#   3. Enrichment testing + plots. comethyl only labels regions; it
#      never tests for enrichment. moduleAnnotatrEnrichment() adds
#      per-module 2x2 Fisher's exact tests against a stated background.
#
# Gene assignment (gene_symbol) is NOT done here. It is done by the
# GREAT functions further down (12a default, --gene_assignment great,
# matching comethyl), or by offline_nearest_gene() with
# --gene_assignment nearest.
# ================================================================

# ------------------------------------------------------------
# Category key: annotatr type suffix -> display label, flag column,
# comethyl-style label (as in comethyl's gene_context/CpG_context
# strings), and family. Single source of truth for everything below.
# ------------------------------------------------------------
annotatr_category_table <- function() {
  data.frame(
    suffix = c("genes_promoters", "genes_1to5kb", "genes_5UTRs", "genes_exons",
               "genes_introns", "genes_3UTRs", "genes_intergenic", "enhancers_fantom",
               "cpg_islands", "cpg_shores", "cpg_shelves", "cpg_inter"),
    label = c("Promoter", "1-5kb upstream", "5' UTR", "Exon",
              "Intron", "3' UTR", "Intergenic", "Enhancer",
              "Island", "Shore", "Shelf", "Open sea"),
    column = c("genic_promoter", "genic_1to5kb", "genic_5UTR", "genic_exon",
               "genic_intron", "genic_3UTR", "genic_intergenic", "genic_enhancer",
               "cpg_island", "cpg_shore", "cpg_shelf", "cpg_opensea"),
    comethyl_label = c("promoter", "1to5kb", "5UTR", "exon",
                       "intron", "3UTR", "intergenic", "enhancer",
                       "island", "shore", "shelf", "open sea"),
    family = c(rep("genic", 8), rep("cpg", 4)),
    stringsAsFactors = FALSE
  )
}

# Default precedence for the single "primary" label used by stacked
# distribution plots (annotatr itself reports every overlap). Enhancer
# is deliberately NOT in the genic precedence: enhancers overlap
# introns/intergenic DNA, so they are an orthogonal axis and are shown
# as an overlap percentage instead.
annotatr_default_priority <- function() {
  list(
    genic = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "1-5kb upstream", "Intergenic"),
    cpg   = c("Island", "Shore", "Shelf", "Open sea")
  )
}

# ------------------------------------------------------------
# Chromosome lengths for the genome, taken from the TxDb package so it
# works offline. Needed so gaps() extends open-sea / intergenic ranges
# to the chromosome ends (without lengths, regions beyond the last
# feature on a chromosome get no CpG label).
# ------------------------------------------------------------
annotatr_genome_seqinfo <- function(genome = "hg38", txdb_pkg = NULL, verbose = TRUE) {
  if (is.null(txdb_pkg)) txdb_pkg <- resolve_annotation_packages(genome)$txdb
  if (requireNamespace(txdb_pkg, quietly = TRUE)) {
    txdb <- get(txdb_pkg, envir = asNamespace(txdb_pkg))
    return(GenomeInfoDb::seqinfo(txdb))
  }
  si <- tryCatch(GenomeInfoDb::Seqinfo(genome = genome), error = function(e) NULL)
  if (is.null(si)) {
    warning("[annotatr] Could not obtain chromosome lengths for ", genome,
            " (", txdb_pkg, " not installed and no internet). Open-sea ranges will ",
            "stop at the last CpG feature on each chromosome.", call. = FALSE)
  }
  si
}

# ------------------------------------------------------------
# .annotatr_std_mcols
#   Internal: strip mcols from a GRanges and add annotatr's standard
#   columns (id = "<prefix>:<n>", tx_id/gene_id/symbol = NA, type) so
#   locally built features combine with annotatr's own.
# ------------------------------------------------------------
.annotatr_std_mcols <- function(gr, type, id_prefix) {
  gr <- GenomicRanges::granges(gr)
  GenomicRanges::mcols(gr)$id      <- paste0(id_prefix, ":", seq_along(gr))
  GenomicRanges::mcols(gr)$tx_id   <- NA_character_
  GenomicRanges::mcols(gr)$gene_id <- NA_character_
  GenomicRanges::mcols(gr)$symbol  <- NA_character_
  GenomicRanges::mcols(gr)$type    <- type
  gr
}

# ------------------------------------------------------------
# CpG islands/shores/shelves/inter from a LOCAL UCSC cpgIslandExt file,
# using annotatr::build_cpg_annots()'s exact algorithm (annotatr 1.28):
#   shores  = flank(islands, 2000) both sides, reduce, setdiff(islands)
#   shelves = flank(shores, 2000) both sides, reduce,
#             setdiff(islands), setdiff(shores)
#   inter   = gaps(union(islands, shores, shelves)), "*" strand only,
#             chromosomes that carry islands only
# annotatr downloads the same UCSC hg38 cpgIslandExt table, so with
# the same file this gives identical ranges.
# ------------------------------------------------------------
build_cpg_annotations_annotatr_local <- function(cpg_island_file, genome = "hg38",
                                                 seqinfo = NULL, verbose = TRUE) {
  validate_file_exists(cpg_island_file, "cpg_island_file")
  cols <- c("bin", "chrom", "chromStart", "chromEnd", "name", "length",
            "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
  df <- utils::read.table(cpg_island_file, sep = "\t", stringsAsFactors = FALSE,
                          col.names = cols, quote = "", comment.char = "")
  if (verbose) message("[annotatr-local] Loaded ", nrow(df), " CpG islands from ", cpg_island_file)

  # UCSC 0-based half-open -> 1-based closed
  if (!is.null(seqinfo)) {
    keep <- df$chrom %in% GenomeInfoDb::seqnames(seqinfo)
    df <- df[keep, , drop = FALSE]
    islands <- GenomicRanges::GRanges(df$chrom, IRanges::IRanges(df$chromStart + 1L, df$chromEnd),
                                      strand = "*", seqinfo = seqinfo)
  } else {
    islands <- GenomicRanges::GRanges(df$chrom, IRanges::IRanges(df$chromStart + 1L, df$chromEnd),
                                      strand = "*")
  }
  islands <- GenomicRanges::sort(islands)

  shores <- c(GenomicRanges::flank(islands, width = 2000, start = TRUE,  both = FALSE),
              GenomicRanges::flank(islands, width = 2000, start = FALSE, both = FALSE))
  shores <- GenomicRanges::reduce(GenomicRanges::trim(GenomicRanges::sort(shores)))
  shores <- GenomicRanges::setdiff(shores, islands)

  shelves <- c(GenomicRanges::flank(shores, width = 2000, start = TRUE,  both = FALSE),
               GenomicRanges::flank(shores, width = 2000, start = FALSE, both = FALSE))
  shelves <- GenomicRanges::reduce(GenomicRanges::trim(GenomicRanges::sort(shelves)))
  shelves <- GenomicRanges::setdiff(shelves, islands)
  shelves <- GenomicRanges::setdiff(shelves, shores)

  extended <- Reduce(function(x, y) GenomicRanges::union(x, y), list(islands, shores, shelves))
  inter <- GenomicRanges::gaps(extended)
  inter <- inter[as.character(GenomicRanges::strand(inter)) == "*"]
  inter <- inter[as.character(GenomicRanges::seqnames(inter)) %in% GenomeInfoDb::seqlevelsInUse(islands)]

  out <- c(.annotatr_std_mcols(islands, paste0(genome, "_cpg_islands"), "island"),
           .annotatr_std_mcols(shores,  paste0(genome, "_cpg_shores"),  "shore"),
           .annotatr_std_mcols(shelves, paste0(genome, "_cpg_shelves"), "shelf"),
           .annotatr_std_mcols(inter,   paste0(genome, "_cpg_inter"),   "inter"))
  if (verbose) {
    message("[annotatr-local] islands=", length(islands), " shores=", length(shores),
            " shelves=", length(shelves), " inter=", length(inter))
  }
  out
}

# ------------------------------------------------------------
# FANTOM5 enhancers from a LOCAL BED. annotatr's hg38 enhancers are the
# FANTOM5 phase 1+2 permissive enhancer set (hg19 BED,
# human_permissive_enhancers_phase_1_and_2.bed.gz) lifted to hg38 with
# the UCSC hg19ToHg38 chain (AnnotationHub AH14150). To reproduce
# offline: download that BED and hg19ToHg38.over.chain.gz once, pass
# both. If the BED is already in the target genome, omit the chain.
# ------------------------------------------------------------
build_enhancer_annotations_local <- function(enhancer_file, genome = "hg38",
                                             chain_file = NULL, seqinfo = NULL, verbose = TRUE) {
  validate_file_exists(enhancer_file, "enhancer_file")
  enh <- rtracklayer::import.bed(enhancer_file)
  if (verbose) message("[annotatr-local] Loaded ", length(enh), " enhancers from ", enhancer_file)
  if (!is.null(chain_file)) {
    validate_file_exists(chain_file, "enhancer_chain_file")
    chain_path <- chain_file
    if (grepl("\\.gz$", chain_file)) {
      chain_path <- tempfile(fileext = ".chain")
      con_in <- gzfile(chain_file, "rb"); on.exit(close(con_in), add = TRUE)
      writeLines(readLines(con_in), chain_path)
    }
    chain <- rtracklayer::import.chain(chain_path)
    enh <- unlist(rtracklayer::liftOver(enh, chain))
    if (verbose) message("[annotatr-local] ", length(enh), " enhancer ranges after liftOver")
  }
  enh <- GenomicRanges::sort(enh)
  if (!is.null(seqinfo)) {
    enh <- GenomeInfoDb::keepSeqlevels(enh, intersect(GenomeInfoDb::seqlevels(enh),
                                                      GenomeInfoDb::seqlevels(seqinfo)),
                                       pruning.mode = "coarse")
  }
  .annotatr_std_mcols(enh, paste0(genome, "_enhancers_fantom"), "enhancer")
}

# ------------------------------------------------------------
# Build (or load) the full annotatr annotation set used by
# comethyl::annotateModule(). Returns a GRanges with a `type` mcol and a
# "provenance" attribute recording where each component came from.
#
# Sources, in order of preference:
#   genes/intergenic: annotatr::build_annotations() (offline; TxDb+org.db)
#   CpGs:      cpg_island_file if given, else annotatr (needs internet)
#   enhancers: enhancer_file if given, else annotatr (needs internet)
# If cache_file exists (and rebuild = FALSE) it is loaded instead.
# ------------------------------------------------------------
build_annotatr_annotations <- function(genome = "hg38",
                                       include_enhancers = TRUE,
                                       cache_file = NULL,
                                       rebuild = FALSE,
                                       cpg_island_file = NULL,
                                       enhancer_file = NULL,
                                       enhancer_chain_file = NULL,
                                       verbose = TRUE) {
  need_types <- paste0(genome, "_", annotatr_category_table()$suffix)
  if (!isTRUE(include_enhancers)) need_types <- setdiff(need_types, paste0(genome, "_enhancers_fantom"))

  if (!is.null(cache_file) && file.exists(cache_file) && !isTRUE(rebuild)) {
    if (verbose) message("[annotatr] Loading cached annotations: ", cache_file)
    annots <- readRDS(cache_file)
    missing_types <- setdiff(need_types, unique(annots$type))
    if (length(missing_types) > 0) {
      stop("Cached annotations in ", cache_file, " are missing: ",
           paste(missing_types, collapse = ", "),
           ". Rebuild the cache (--rebuild_annotatr_cache TRUE), or set ",
           "--include_enhancers FALSE if only enhancers are missing.", call. = FALSE)
    }
    return(annots)
  }

  if (!requireNamespace("annotatr", quietly = TRUE)) {
    stop("annotatr is not installed in this environment. Add bioconductor-annotatr to the pixi environment.",
         call. = FALSE)
  }
  seqinfo <- annotatr_genome_seqinfo(genome, verbose = verbose)
  provenance <- list(annotatr_version = as.character(utils::packageVersion("annotatr")),
                     genome = genome, built = timestamp_now())

  # Genes + intergenic (offline: TxDb + org.db)
  if (verbose) message("[annotatr] Building ", genome, "_basicgenes + ", genome, "_genes_intergenic")
  genes <- suppressWarnings(suppressMessages(
    annotatr::build_annotations(genome = genome,
                                annotations = paste0(genome, c("_basicgenes", "_genes_intergenic")))
  ))
  provenance$genes <- paste0("annotatr::build_annotations(", genome, "_basicgenes, ",
                             genome, "_genes_intergenic); TxDb ",
                             resolve_annotation_packages(genome)$txdb, " ",
                             as.character(utils::packageVersion(resolve_annotation_packages(genome)$txdb)))

  # CpGs
  if (!is.null(cpg_island_file)) {
    cpgs <- build_cpg_annotations_annotatr_local(cpg_island_file, genome = genome,
                                                 seqinfo = seqinfo, verbose = verbose)
    provenance$cpgs <- paste0("local file (annotatr algorithm): ", cpg_island_file,
                              " md5=", unname(tools::md5sum(cpg_island_file)))
  } else {
    if (verbose) message("[annotatr] Building ", genome, "_cpgs via annotatr (requires internet)")
    cpgs <- suppressWarnings(suppressMessages(
      annotatr::build_annotations(genome = genome, annotations = paste0(genome, "_cpgs"))
    ))
    provenance$cpgs <- paste0("annotatr::build_annotations(", genome, "_cpgs) [UCSC cpgIslandExt download]")
  }

  pieces <- list(genes, cpgs)

  # Enhancers
  if (isTRUE(include_enhancers)) {
    if (!is.null(enhancer_file)) {
      enh <- build_enhancer_annotations_local(enhancer_file, genome = genome,
                                              chain_file = enhancer_chain_file,
                                              seqinfo = seqinfo, verbose = verbose)
      provenance$enhancers <- paste0("local file: ", enhancer_file,
                                     ifelse(is.null(enhancer_chain_file), "",
                                            paste0(" lifted with ", enhancer_chain_file)))
    } else {
      if (verbose) message("[annotatr] Building ", genome, "_enhancers_fantom via annotatr (requires internet)")
      enh <- suppressWarnings(suppressMessages(
        annotatr::build_annotations(genome = genome, annotations = paste0(genome, "_enhancers_fantom"))
      ))
      provenance$enhancers <- paste0("annotatr::build_annotations(", genome,
                                     "_enhancers_fantom) [FANTOM5 phase 1+2 permissive, hg19 lifted via AH14150 for hg38]")
    }
    pieces <- c(pieces, list(enh))
  }

  # Harmonise mcols before combining (annotatr returns id/tx_id/gene_id/symbol/type)
  pieces <- lapply(pieces, function(g) {
    g <- GenomicRanges::granges(g, use.mcols = TRUE)
    keep <- c("id", "tx_id", "gene_id", "symbol", "type")
    for (k in setdiff(keep, colnames(GenomicRanges::mcols(g)))) GenomicRanges::mcols(g)[[k]] <- NA_character_
    GenomicRanges::mcols(g) <- GenomicRanges::mcols(g)[, keep]
    for (k in keep) GenomicRanges::mcols(g)[[k]] <- as.character(GenomicRanges::mcols(g)[[k]])
    g
  })
  annots <- do.call(c, unname(pieces))  # c() merges compatible seqinfo

  # As comethyl::annotateModule(): standard chromosomes only
  annots <- GenomeInfoDb::keepStandardChromosomes(annots, pruning.mode = "coarse")

  missing_types <- setdiff(need_types, unique(annots$type))
  if (length(missing_types) > 0) {
    stop("annotatr build is missing annotation types: ", paste(missing_types, collapse = ", "),
         call. = FALSE)
  }
  attr(annots, "provenance") <- provenance

  if (!is.null(cache_file)) {
    safe_dir_create(dirname(cache_file))
    saveRDS(annots, cache_file)
    if (verbose) message("[annotatr] Saved annotation cache: ", cache_file)
  }
  annots
}

# ------------------------------------------------------------
# Label every region with every annotatr feature it overlaps.
#
# Equivalent to annotatr::annotate_regions(regions, annotations,
# minoverlap = 1, ignore.strand = TRUE) followed by comethyl's
# aggregation into comma-separated strings, but computed with
# overlapsAny() per annotation type. annotate_regions() returns one
# row per region x overlapping transcript feature, which for ~465k
# regions is tens of millions of rows; overlapsAny() gives the same
# yes/no answer per region and category without that blow-up.
#
# Adds:
#   one logical flag column per category (genic_*, cpg_*)
#   gene_context / CpG_context -- comethyl-style strings,
#       e.g. "exon, intron" / "shore"
#   gene_context_primary / CpG_context_primary -- single label by
#       precedence (see annotatr_default_priority())
# ------------------------------------------------------------
annotate_regions_annotatr <- function(regions_df, annotations, genome = "hg38",
                                      priority = annotatr_default_priority(),
                                      verbose = TRUE) {
  validate_regions_df(regions_df, "regions_df", req_cols = c("RegionID", "chr", "start", "end"))
  cat_tab <- annotatr_category_table()
  types_present <- unique(annotations$type)
  cat_tab <- cat_tab[paste0(genome, "_", cat_tab$suffix) %in% types_present, , drop = FALSE]

  gr <- GenomicRanges::GRanges(regions_df$chr, IRanges::IRanges(regions_df$start, regions_df$end))

  out <- regions_df
  for (i in seq_len(nrow(cat_tab))) {
    ty <- paste0(genome, "_", cat_tab$suffix[i])
    out[[cat_tab$column[i]]] <- IRanges::overlapsAny(
      gr, annotations[annotations$type == ty], minoverlap = 1L, ignore.strand = TRUE
    )
    if (verbose) message("[annotatr] ", sprintf("%-15s", cat_tab$label[i]), " ",
                         sum(out[[cat_tab$column[i]]]), " regions")
  }

  collapse_labels <- function(cols, labels) {
    m <- as.matrix(out[, cols, drop = FALSE])
    s <- apply(m, 1, function(x) paste(labels[x], collapse = ", "))
    s[!nzchar(s)] <- NA_character_
    s
  }
  pick_primary <- function(cols, labels, order) {
    res <- rep(NA_character_, nrow(out))
    for (lab in rev(order)) {           # last assignment wins -> highest priority last
      j <- match(lab, labels)
      if (is.na(j)) next
      res[out[[cols[j]]]] <- lab
    }
    res
  }

  g <- cat_tab[cat_tab$family == "genic", ]
  c_ <- cat_tab[cat_tab$family == "cpg", ]
  out$gene_context <- collapse_labels(g$column, g$comethyl_label)
  out$CpG_context  <- collapse_labels(c_$column, c_$comethyl_label)
  out$gene_context_primary <- pick_primary(g$column, g$label, priority$genic)
  out$CpG_context_primary  <- pick_primary(c_$column, c_$label, priority$cpg)

  n_no_cpg <- sum(is.na(out$CpG_context))
  if (n_no_cpg > 0 && verbose) {
    message("[annotatr] NOTE: ", n_no_cpg, " regions overlap no CpG category ",
            "(non-standard chromosome, or open sea built without chromosome lengths)")
  }
  out
}

# ------------------------------------------------------------
# Per-module enrichment of annotatr categories: 2x2 Fisher's exact
# test of (in module vs not in module) x (overlaps category vs not),
# within a stated background.
#
# background = "all"      -> every region passed into network
#                            construction, grey included (default;
#                            matches comethyl::enrichModule()'s GREAT
#                            background and DMRichR's "all testable
#                            regions"). "Not in module" = all other
#                            regions, grey included.
# background = "assigned" -> only module-assigned (non-grey) regions;
#                            each module vs the rest of the network.
#
# include_pooled adds one extra test, "All modules" = every non-grey
# region vs grey (only meaningful with background = "all").
#
# adjust_width = TRUE adds a sensitivity estimate from a logistic
# regression  overlap ~ in_module + log(width), because wider regions
# are more likely to overlap any feature (strongly so for CpG islands).
#
# FDR = Benjamini-Hochberg across all tests in the call (one family).
# ------------------------------------------------------------
moduleAnnotatrEnrichment <- function(annotated_regions,
                                     family = c("genic", "cpg"),
                                     background = c("all", "assigned"),
                                     exclude_modules = "grey",
                                     include_pooled = TRUE,
                                     pooled_label = "All modules",
                                     adjust_width = FALSE,
                                     module_col = "module",
                                     verbose = TRUE) {
  family <- match.arg(family)
  background <- match.arg(background)
  cat_tab <- annotatr_category_table()
  cat_tab <- cat_tab[cat_tab$family == family & cat_tab$column %in% colnames(annotated_regions), , drop = FALSE]
  if (nrow(cat_tab) == 0) {
    stop("moduleAnnotatrEnrichment(): no ", family, " flag columns found -- run annotate_regions_annotatr() first.",
         call. = FALSE)
  }

  df <- annotated_regions
  mod_vec <- as.character(df[[module_col]])
  if (background == "assigned") {
    keep <- !(mod_vec %in% exclude_modules)
    df <- df[keep, , drop = FALSE]; mod_vec <- mod_vec[keep]
    if (isTRUE(include_pooled) && verbose) {
      message("[enrich] background = 'assigned': skipping pooled '", pooled_label,
              "' test (it would be the whole background)")
    }
    include_pooled <- FALSE
  }
  modules <- sort(unique(mod_vec[!(mod_vec %in% exclude_modules)]))

  logw <- NULL
  if (isTRUE(adjust_width)) {
    w <- if ("width" %in% colnames(df)) df$width else (df$end - df$start + 1)
    logw <- log(pmax(as.numeric(w), 1))
  }

  test_sets <- lapply(modules, function(m) mod_vec == m); names(test_sets) <- modules
  if (isTRUE(include_pooled)) {
    test_sets <- c(stats::setNames(list(!(mod_vec %in% exclude_modules)), pooled_label), test_sets)
  }

  if (verbose) message("[enrich] ", family, ": ", length(test_sets), " sets x ", nrow(cat_tab),
                       " categories; background = ", background, " (n = ", nrow(df), ")")

  width_adjusted <- function(y, x) {
    na_out <- c(OR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
    if (sum(y & x) == 0 || sum(!y & x) == 0) return(na_out)  # complete separation in module
    X <- cbind(1, as.numeric(x), logw)
    fit <- tryCatch(suppressWarnings(stats::glm.fit(X, as.numeric(y), family = stats::binomial())),
                    error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged)) return(na_out)
    pr <- fit$fitted.values
    V <- tryCatch(solve(crossprod(X * sqrt(pr * (1 - pr)))), error = function(e) NULL)
    if (is.null(V)) return(na_out)
    b <- fit$coefficients[2]; se <- sqrt(V[2, 2])
    c(OR = exp(b), lo = exp(b - 1.959964 * se), hi = exp(b + 1.959964 * se),
      p = 2 * stats::pnorm(-abs(b / se)))
  }

  rows <- vector("list", length(test_sets) * nrow(cat_tab)); k <- 0
  for (set_name in names(test_sets)) {
    in_set <- test_sets[[set_name]]
    for (i in seq_len(nrow(cat_tab))) {
      has <- as.logical(df[[cat_tab$column[i]]]); has[is.na(has)] <- FALSE
      a <- sum(in_set & has);  b <- sum(in_set & !has)
      c <- sum(!in_set & has); d <- sum(!in_set & !has)
      ft <- stats::fisher.test(matrix(c(a, c, b, d), nrow = 2))
      row <- data.frame(
        module = set_name, category = cat_tab$label[i],
        n_module = a + b, n_in_module = a,
        pct_in_module = 100 * a / (a + b),
        pct_rest = 100 * c / (c + d),
        pct_background = 100 * (a + c) / (a + b + c + d),
        OR = unname(ft$estimate), CIlower = ft$conf.int[1], CIupper = ft$conf.int[2],
        p = ft$p.value, stringsAsFactors = FALSE
      )
      if (isTRUE(adjust_width)) {
        wa <- width_adjusted(has, in_set)
        row$OR_width_adj <- wa[["OR"]]; row$CIlower_width_adj <- wa[["lo"]]
        row$CIupper_width_adj <- wa[["hi"]]; row$p_width_adj <- wa[["p"]]
      }
      k <- k + 1; rows[[k]] <- row
    }
  }
  out <- do.call(rbind, rows[seq_len(k)])
  out$fdr <- stats::p.adjust(out$p, method = "BH")
  if (isTRUE(adjust_width)) out$fdr_width_adj <- stats::p.adjust(out$p_width_adj, method = "BH")
  out$background <- background
  out$n_background <- nrow(df)
  out$family <- family
  front <- c("family", "module", "category", "n_module", "n_in_module", "pct_in_module",
             "pct_rest", "pct_background", "OR", "CIlower", "CIupper", "p", "fdr")
  out[, c(front, setdiff(colnames(out), front))]
}

# ------------------------------------------------------------
# Distribution table: for each group (Background, All modules, each
# module) and category, the % of regions overlapping that category
# (annotatr multi-label; columns do not sum to 100) and the % with that
# category as their primary label (sums to 100 within family).
# ------------------------------------------------------------
summarizeAnnotationDistribution <- function(annotated_regions,
                                            exclude_modules = "grey",
                                            background_label = "Background",
                                            pooled_label = "All modules",
                                            module_col = "module") {
  cat_tab <- annotatr_category_table()
  cat_tab <- cat_tab[cat_tab$column %in% colnames(annotated_regions), , drop = FALSE]
  mod_vec <- as.character(annotated_regions[[module_col]])
  groups <- c(list(rep(TRUE, length(mod_vec)), !(mod_vec %in% exclude_modules)),
              lapply(sort(unique(mod_vec[!(mod_vec %in% exclude_modules)])), function(m) mod_vec == m))
  names(groups) <- c(background_label, pooled_label,
                     sort(unique(mod_vec[!(mod_vec %in% exclude_modules)])))
  prim <- list(genic = annotated_regions$gene_context_primary, cpg = annotated_regions$CpG_context_primary)

  rows <- list()
  for (gname in names(groups)) {
    idx <- groups[[gname]]; n <- sum(idx)
    for (i in seq_len(nrow(cat_tab))) {
      fam <- cat_tab$family[i]; lab <- cat_tab$label[i]
      n_ov <- sum(annotated_regions[[cat_tab$column[i]]][idx], na.rm = TRUE)
      n_pr <- sum(prim[[fam]][idx] == lab, na.rm = TRUE)
      rows[[length(rows) + 1]] <- data.frame(
        group = gname, family = fam, category = lab, n_group = n,
        n_overlap = n_ov, pct_overlap = 100 * n_ov / n,
        n_primary = n_pr, pct_primary = 100 * n_pr / n,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$pct_primary[out$category == "Enhancer"] <- NA  # enhancer not in genic precedence
  out$n_primary[out$category == "Enhancer"] <- NA
  out
}

# ------------------------------------------------------------
# Plot helpers
# ------------------------------------------------------------
# ------------------------------------------------------------
# annot_parse_color_overrides
#   Parse "yellow=darkgoldenrod3,white=grey80" into a named vector
#   (names = module, values = colour). Empty/NULL -> character(0).
# ------------------------------------------------------------
annot_parse_color_overrides <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  pairs <- strsplit(strsplit(x, ",")[[1]], "=")
  bad <- vapply(pairs, length, 1L) != 2
  if (any(bad)) stop("--color_overrides must look like 'yellow=darkgoldenrod3,white=grey80'", call. = FALSE)
  stats::setNames(trimws(vapply(pairs, `[`, "", 2)), trimws(vapply(pairs, `[`, "", 1)))
}

# ------------------------------------------------------------
# annot_group_colors
#   Bar colour per group: override if given; grey35 for Background,
#   grey65 for All modules; the module name itself if it is a valid R
#   colour; else grey50. Returns a named vector.
# ------------------------------------------------------------
annot_group_colors <- function(groups, color_overrides = character(0),
                               background_label = "Background", pooled_label = "All modules") {
  is_col <- function(x) !inherits(tryCatch(grDevices::col2rgb(x), error = function(e) e), "error")
  cols <- vapply(groups, function(g) {
    if (g %in% names(color_overrides)) return(color_overrides[[g]])
    if (g == background_label) return("grey35")
    if (g == pooled_label) return("grey65")
    if (is_col(g)) return(g)
    "grey50"
  }, "")
  stats::setNames(cols, groups)
}

# ------------------------------------------------------------
# annot_category_palette
#   Fixed colour palette for annotatr categories in one family
#   ("genic" includes Enhancer; "cpg" = Island/Shore/Shelf/Open sea).
# ------------------------------------------------------------
annot_category_palette <- function(family) {
  if (family == "genic") {
    c("Promoter" = "#D55E00", "5' UTR" = "#E69F00", "3' UTR" = "#CC79A7", "Exon" = "#009E73",
      "Intron" = "#56B4E9", "1-5kb upstream" = "#F0E442", "Intergenic" = "#BDBDBD", "Enhancer" = "#0072B2")
  } else {
    c("Island" = "#00441B", "Shore" = "#238B45", "Shelf" = "#74C476", "Open sea" = "#E5F5E0")
  }
}

# ------------------------------------------------------------
# save_annotation_plot
#   Save a ggplot as <file_stem>.pdf (cairo_pdf when available) and,
#   if png = TRUE, <file_stem>.png at `dpi`. Returns the paths invisibly.
# ------------------------------------------------------------
save_annotation_plot <- function(p, file_stem, width, height, dpi = 600, verbose = TRUE, png = TRUE) {
  safe_dir_create(dirname(file_stem))
  pdf_file <- paste0(file_stem, ".pdf"); png_file <- paste0(file_stem, ".png")
  dev <- if (capabilities("cairo")) grDevices::cairo_pdf else {
    warning("cairo not available; falling back to pdf() (fonts may not embed as expected)", call. = FALSE)
    grDevices::pdf
  }
  ggplot2::ggsave(pdf_file, p, width = width, height = height, units = "in", device = dev, limitsize = FALSE)
  if (isTRUE(png)) ggplot2::ggsave(png_file, p, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)
  if (verbose) message("[plot] Saved: ", pdf_file, if (isTRUE(png)) " (+ .png)" else "")
  invisible(c(pdf_file, if (isTRUE(png)) png_file))
}

# ------------------------------------------------------------
# .annot_group_levels
#   Internal: x-axis order for distribution plots -- Background and
#   All modules first (if requested), then `groups` in the given order,
#   or all modules largest-first when groups = NULL. Warns on unknown
#   modules.
# ------------------------------------------------------------
.annot_group_levels <- function(dist_df, groups, background_label, pooled_label,
                                include_background, include_pooled) {
  mods <- setdiff(unique(dist_df$group), c(background_label, pooled_label))
  if (!is.null(groups)) {
    missing <- setdiff(groups, mods)
    if (length(missing) > 0) warning("Modules not found and skipped: ", paste(missing, collapse = ", "), call. = FALSE)
    mods <- intersect(groups, mods)
  } else {
    sizes <- dist_df[!duplicated(dist_df$group), c("group", "n_group")]
    sizes <- sizes[sizes$group %in% mods, ]
    mods <- sizes$group[order(-sizes$n_group, sizes$group)]
  }
  c(if (include_background) background_label, if (include_pooled) pooled_label, mods)
}

# ------------------------------------------------------------
# plotAnnotationDistribution()
#   mode = "stacked": 100% stacked bars of the PRIMARY label per group
#          (one label per region by precedence; enhancer excluded).
#   mode = "overlap": % of regions overlapping each category (annotatr
#          multi-label), one facet per category, bars coloured by
#          module, dashed line = background %.
# groups = NULL -> all modules (largest first); otherwise that subset,
# in that order. Background and "All modules" bars are prepended.
# ------------------------------------------------------------
plotAnnotationDistribution <- function(dist_df,
                                       family = c("genic", "cpg"),
                                       mode = c("stacked", "overlap"),
                                       groups = NULL,
                                       include_background = TRUE,
                                       include_pooled = TRUE,
                                       background_label = "Background",
                                       pooled_label = "All modules",
                                       color_overrides = character(0),
                                       show_n = TRUE,
                                       ylim_max = "auto",
                                       title = NULL,
                                       file_stem = NULL,
                                       width = NULL, height = NULL,
                                       base_size = 10,
                                       verbose = TRUE) {
  # ylim_max (stacked mode only): "auto" zooms the CpG stacked plot to the
  # non-open-sea part (open sea is ~95% of regions and hides everything
  # else); NULL = full 0-100%; a number = explicit upper limit in %.
  family <- match.arg(family); mode <- match.arg(mode)
  d <- dist_df[dist_df$family == family, , drop = FALSE]
  lv <- .annot_group_levels(d, groups, background_label, pooled_label, include_background, include_pooled)
  d <- d[d$group %in% lv, , drop = FALSE]

  n_lab <- d[!duplicated(d$group), c("group", "n_group")]
  xlab_map <- if (show_n) {
    stats::setNames(paste0(n_lab$group, "\n(n=", format(n_lab$n_group, big.mark = ",", trim = TRUE), ")"), n_lab$group)
  } else stats::setNames(n_lab$group, n_lab$group)
  d$group <- factor(d$group, levels = lv)
  pal <- annot_category_palette(family)
  fam_title <- if (family == "genic") "Genic annotation" else "CpG context"

  if (mode == "stacked") {
    order_lv <- annotatr_default_priority()[[family]]
    d <- d[d$category %in% order_lv, , drop = FALSE]
    d$category <- factor(d$category, levels = rev(order_lv))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = pct_primary, fill = category)) +
      ggplot2::geom_col(width = 0.8, color = "white", linewidth = 0.15) +
      ggplot2::scale_fill_manual(values = pal, breaks = order_lv, name = NULL) +
      ggplot2::scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
      ggplot2::scale_x_discrete(labels = xlab_map) +
      ggplot2::labs(x = NULL, y = "Regions (%)",
                    title = if (is.null(title)) paste0(fam_title, " (primary label per region)") else title) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = base_size - 2),
                     legend.position = "right")
    ymax <- ylim_max
    if (identical(ymax, "auto")) {
      ymax <- NULL
      if (family == "cpg") {
        top <- order_lv[length(order_lv)]  # "Open sea"
        rest <- tapply(d$pct_primary[d$category != top], d$group[d$category != top], sum, na.rm = TRUE)
        ymax <- min(100, 5 * ceiling(1.15 * max(rest, na.rm = TRUE) / 5))
      }
    }
    if (!is.null(ymax) && ymax < 100) {
      p <- p + ggplot2::coord_cartesian(ylim = c(0, ymax)) +
        ggplot2::labs(subtitle = paste0("Axis cut at ", ymax, "%; rest = ",
                                        tolower(order_lv[length(order_lv)])))
    }
    if (is.null(width))  width  <- max(5, 1.6 + 0.32 * length(lv))
    if (is.null(height)) height <- 4.5
  } else {
    cats <- names(pal)[names(pal) %in% unique(d$category)]
    d$category <- factor(d$category, levels = cats)
    bg <- d[d$group == background_label, c("category", "pct_overlap")]
    if (nrow(bg) == 0) {
      bg <- dist_df[dist_df$family == family & dist_df$group == background_label, c("category", "pct_overlap")]
      bg$category <- factor(bg$category, levels = cats)
    }
    gcols <- annot_group_colors(lv, color_overrides, background_label, pooled_label)
    p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = pct_overlap, fill = group)) +
      ggplot2::geom_col(width = 0.8, color = "grey20", linewidth = 0.15) +
      ggplot2::geom_hline(data = bg, ggplot2::aes(yintercept = pct_overlap),
                          linetype = "dashed", linewidth = 0.3, color = "grey20") +
      ggplot2::facet_wrap(~ category, scales = "free_y", nrow = if (family == "genic") 2 else 1) +
      ggplot2::scale_fill_manual(values = gcols, guide = "none") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08)),
                                  labels = function(x) paste0(x, "%")) +
      ggplot2::scale_x_discrete(labels = if (length(lv) <= 12) xlab_map else ggplot2::waiver()) +
      ggplot2::labs(x = NULL, y = "Regions overlapping category (%)",
                    title = if (is.null(title)) paste0(fam_title, " (% overlapping; dashed = background)") else title) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = base_size - 3),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face = "bold"))
    ncol_f <- if (family == "genic") 4 else 4
    if (is.null(width))  width  <- max(6, min(26, ncol_f * (0.9 + 0.17 * length(lv))))
    if (is.null(height)) height <- if (family == "genic") 6.5 else 3.6
  }
  if (!is.null(file_stem)) save_annotation_plot(p, file_stem, width, height, verbose = verbose)
  invisible(p)
}

# ------------------------------------------------------------
# Enrichment heatmap: module x category, fill = log2(OR), "*" = FDR <
# p_cut. Grey tiles = OR not estimable (0 or Inf: zero counts in a
# small module). log2OR_cap squishes extreme values so one very large
# OR does not wash out the colour scale.
# ------------------------------------------------------------
plotAnnotatrEnrichmentHeatmap <- function(enrich_df,
                                          groups = NULL,
                                          p_cut = 0.05,
                                          log2OR_cap = 3,
                                          pooled_label = "All modules",
                                          use_width_adjusted = FALSE,
                                          title = NULL,
                                          file_stem = NULL,
                                          width = NULL, height = NULL,
                                          base_size = 10,
                                          verbose = TRUE) {
  d <- enrich_df
  if (isTRUE(use_width_adjusted)) {
    if (!"OR_width_adj" %in% colnames(d)) stop("No width-adjusted columns; run with adjust_width = TRUE", call. = FALSE)
    d$OR <- d$OR_width_adj; d$fdr <- d$fdr_width_adj
  }
  fam <- unique(d$family)[1]
  mods <- setdiff(unique(d$module), pooled_label)
  if (!is.null(groups)) mods <- intersect(groups, mods) else {
    sz <- d[!duplicated(d$module) & d$module %in% mods, c("module", "n_module")]
    mods <- sz$module[order(-sz$n_module, sz$module)]
  }
  lv <- c(if (pooled_label %in% d$module) pooled_label, mods)
  d <- d[d$module %in% lv, , drop = FALSE]
  d$log2OR <- log2(d$OR); d$log2OR[!is.finite(d$log2OR)] <- NA
  d$sig <- !is.na(d$fdr) & d$fdr < p_cut & !is.na(d$log2OR)
  lim <- if (is.null(log2OR_cap)) max(abs(d$log2OR), na.rm = TRUE) else log2OR_cap
  if (!is.finite(lim) || lim == 0) lim <- 1
  cats <- names(annot_category_palette(fam))
  cats <- cats[cats %in% unique(d$category)]
  d$module <- factor(d$module, levels = lv)
  d$category <- factor(d$category, levels = rev(cats))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = module, y = category, fill = log2OR)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(sig, "*", "")), size = 3.5, vjust = 0.75) +
    ggplot2::scale_fill_gradientn(
      name = if (is.null(log2OR_cap)) "log2(OR)" else paste0("log2(OR)\ncapped at +/-", lim),
      colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"),
      limits = c(-lim, lim), oob = scales::squish, na.value = "grey85") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = if (is.null(title)) paste0(if (fam == "genic") "Genic" else "CpG",
                                                     " enrichment by module",
                                                     if (isTRUE(use_width_adjusted)) " (width-adjusted)" else "") else title,
                  caption = paste0("* FDR < ", p_cut, "; grey = OR not estimable; background = ",
                                   unique(d$background)[1], " regions")) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = base_size - 2),
                   panel.grid = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(size = base_size - 3, colour = "grey30"))
  if (is.null(width))  width  <- max(5, min(24, 2.6 + 0.2 * length(lv)))
  if (is.null(height)) height <- max(2.8, 1.4 + 0.32 * length(cats))
  if (!is.null(file_stem)) save_annotation_plot(p, file_stem, width, height, verbose = verbose)
  invisible(p)
}

# ------------------------------------------------------------
# Forest plot for a few modules of interest: log2(OR) with 95% CI per
# category, one panel per module. Filled = FDR < p_cut. Non-estimable
# ORs (0/Inf) are dropped and listed in the caption.
# ------------------------------------------------------------
plotAnnotatrEnrichmentForest <- function(enrich_df,
                                         groups,
                                         p_cut = 0.05,
                                         log2OR_cap = 4,
                                         pooled_label = "All modules",
                                         include_pooled = TRUE,
                                         color_overrides = character(0),
                                         use_width_adjusted = FALSE,
                                         title = NULL,
                                         file_stem = NULL,
                                         width = NULL, height = NULL,
                                         base_size = 10,
                                         verbose = TRUE) {
  d <- enrich_df
  if (isTRUE(use_width_adjusted)) {
    d$OR <- d$OR_width_adj; d$CIlower <- d$CIlower_width_adj
    d$CIupper <- d$CIupper_width_adj; d$fdr <- d$fdr_width_adj
  }
  fam <- unique(d$family)[1]
  lv <- c(if (include_pooled && pooled_label %in% d$module) pooled_label, intersect(groups, unique(d$module)))
  d <- d[d$module %in% lv, , drop = FALSE]
  dropped <- d[!is.finite(log2(d$OR)) | !is.finite(log2(d$CIlower)) | !is.finite(log2(d$CIupper)), ]
  d <- d[is.finite(log2(d$OR)), , drop = FALSE]
  d$l2 <- log2(d$OR)
  d$lo <- pmax(log2(pmax(d$CIlower, .Machine$double.xmin)), -log2OR_cap)
  d$hi <- pmin(log2(d$CIupper), log2OR_cap)
  d$l2 <- pmin(pmax(d$l2, -log2OR_cap), log2OR_cap)
  d$sig <- !is.na(d$fdr) & d$fdr < p_cut
  cats <- names(annot_category_palette(fam)); cats <- cats[cats %in% unique(d$category)]
  d$category <- factor(d$category, levels = rev(cats))
  d$module <- factor(d$module, levels = lv)
  gcols <- annot_group_colors(lv, color_overrides, pooled_label = pooled_label)

  cap_txt <- paste0("Filled = FDR < ", p_cut, "; bars = 95% CI, truncated at +/-", log2OR_cap,
                    " log2 units; background = ", unique(enrich_df$background)[1], " regions")
  if (nrow(dropped) > 0) {
    cap_txt <- paste0(cap_txt, "\nNot estimable (zero counts): ",
                      paste(unique(paste0(dropped$module, "/", dropped$category)), collapse = ", "))
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x = l2, y = category)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(fill = module, shape = sig), size = 2.4, stroke = 0.5) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), guide = "none") +
    ggplot2::scale_fill_manual(values = gcols, guide = "none") +
    ggplot2::facet_wrap(~ module, nrow = 1) +
    ggplot2::coord_cartesian(xlim = c(-log2OR_cap, log2OR_cap)) +
    ggplot2::labs(x = expression(log[2]~"(odds ratio)"), y = NULL,
                  title = if (is.null(title)) paste0(if (fam == "genic") "Genic" else "CpG",
                                                     " enrichment, modules of interest") else title,
                  caption = cap_txt) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.caption = ggplot2::element_text(size = base_size - 3, colour = "grey30", hjust = 0))
  if (is.null(width))  width  <- max(4, 1.4 + 1.9 * length(lv))
  if (is.null(height)) height <- max(2.6, 1.6 + 0.3 * length(cats))
  if (!is.null(file_stem)) save_annotation_plot(p, file_stem, width, height, verbose = verbose)
  invisible(p)
}

# ================================================================
# GREAT REGION -> GENE ASSIGNMENT (offline; matches comethyl::annotateModule())
#
# comethyl::annotateModule() sends regions to the GREAT web server
# (v4.0.4, basalPlusExt: basal domain 5 kb upstream / 1 kb downstream of
# each TSS, extended up to 1 Mb) and gets back up to two genes per
# region plus distance to TSS. The functions below reproduce that with
# rGREAT's LOCAL implementation, so it runs on compute nodes:
#
#   tss_source = "GREAT:hg38" (default) -> GREAT v4's own gene/TSS table,
#       shipped inside rGREAT (extdata/GREATv4.genes.hg38.tsv.gz; the
#       same ~18.8k mostly protein-coding genes the web server uses).
#       Built exactly as rGREAT::getGREATDefaultTSS() does, except that
#       chromosome lengths come from rGREAT's shipped CHR_LEN_DB instead
#       of a live UCSC download.
#   tss_source = "TxDb:hg38" -> rGREAT::extendTSSFromTxDb() on
#       TxDb.Hsapiens.UCSC.hg38.knownGene (protein-coding genes only,
#       one TSS per gene = gene 5' end).
#
# Region -> gene rule = rGREAT::getRegionGeneAssociations(): a region is
# assigned to every gene whose regulatory domain it overlaps (usually 1-2
# genes); distance = signed distance from the nearest region edge to the
# TSS (negative = upstream of the gene, 0 = region covers the TSS).
#
# Gene IDs are Entrez throughout (rGREAT's convention), so gene sets for
# enrichment must be Entrez too.
# ================================================================

great_parse_tss_source <- function(tss_source, genome = "hg38") {
  x <- trimws(tss_source)
  if (grepl("^great(:|$)", x, ignore.case = TRUE)) {
    g <- sub("^great:?", "", x, ignore.case = TRUE); if (!nzchar(g)) g <- genome
    return(list(category = "GREAT", genome = g, label = paste0("GREAT:", g)))
  }
  if (grepl("^txdb(:|$)", x, ignore.case = TRUE)) {
    g <- sub("^txdb:?", "", x, ignore.case = TRUE); if (!nzchar(g)) g <- genome
    return(list(category = "TxDb", genome = g, label = paste0("TxDb:", g),
                txdb_pkg = resolve_annotation_packages(g)$txdb))
  }
  stop("--tss_source must be 'GREAT:<genome>' (default) or 'TxDb:<genome>'; got '", tss_source, "'",
       call. = FALSE)
}

# Ensembl/symbol -> Entrez for the GREAT v4 table, identical to the logic in
# rGREAT::getGREATDefaultTSS() (org.Hs.eg.db ENSEMBL2EG first, SYMBOL2EG as
# fallback; comma-separated Ensembl lists take the first mapped ID).
.great_map_table_to_entrez <- function(tb, genome = "hg38") {
  orgdb <- if (genome %in% c("hg19", "hg38")) "org.Hs.eg.db" else "org.Mm.eg.db"
  if (!requireNamespace(orgdb, quietly = TRUE)) stop(orgdb, " is required for GREAT TSS mapping.", call. = FALSE)
  pre <- if (orgdb == "org.Hs.eg.db") "org.Hs.eg" else "org.Mm.eg"
  map1 <- unlist(as.list(get(paste0(pre, "ENSEMBL2EG"), envir = asNamespace(orgdb))))
  map2 <- unlist(as.list(get(paste0(pre, "SYMBOL2EG"), envir = asNamespace(orgdb))))
  id1 <- map1[tb[, 1]]
  l <- nchar(tb[, 1]) > 15
  id1[l] <- vapply(strsplit(tb[l, 1], ","), function(x) {
    x2 <- map1[x]; x2 <- x2[!is.na(x2)]; if (length(x2) == 0) NA_character_ else x2[[1]]
  }, "")
  l2 <- grep("^ENS.*\\d+$", tb[, 5])
  id2 <- map2[tb[, 5]]
  id2[l2] <- map1[tb[l2, 5]]
  unname(ifelse(is.na(id1), id2, id1))
}

# Entrez -> symbol / description / Ensembl from org.Hs.eg.db (offline,
# versioned; replaces comethyl's live biomaRt lookup).
gene_info_from_orgdb <- function(entrez, genome = "hg38") {
  entrez <- unique(as.character(entrez[!is.na(entrez) & nzchar(entrez)]))
  empty <- data.frame(gene_entrezID = character(0), symbol_orgdb = character(0),
                      gene_description = character(0), gene_ensemblID = character(0))
  if (length(entrez) == 0) return(empty)
  orgdb_pkg <- resolve_annotation_packages(genome)$annoDb
  if (!requireNamespace(orgdb_pkg, quietly = TRUE)) stop(orgdb_pkg, " is required.", call. = FALSE)
  orgdb <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))
  m <- suppressMessages(AnnotationDbi::select(orgdb, keys = entrez, keytype = "ENTREZID",
                                              columns = c("SYMBOL", "GENENAME", "ENSEMBL")))
  agg <- function(v) paste(unique(v[!is.na(v)]), collapse = ",")
  out <- do.call(rbind, lapply(split(m, m$ENTREZID), function(d) data.frame(
    gene_entrezID = d$ENTREZID[1], symbol_orgdb = agg(d$SYMBOL),
    gene_description = agg(d$GENENAME), gene_ensemblID = agg(d$ENSEMBL), stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------
# Build the GREAT regulatory domains (rGREAT "extended_tss").
# Returns a GRanges with gene_id (Entrez), tss_position, tss_strand and
# attr "symbol_map" (Entrez -> symbol) and "great_params".
# ------------------------------------------------------------
great_build_extended_tss <- function(tss_source = "GREAT:hg38", genome = "hg38",
                                     mode = c("basalPlusExt", "twoClosest", "oneClosest"),
                                     basal_upstream = 5000, basal_downstream = 1000,
                                     extension = 1e6, verbose = TRUE) {
  mode <- match.arg(mode)
  if (!requireNamespace("rGREAT", quietly = TRUE)) stop("rGREAT is not installed.", call. = FALSE)
  src <- great_parse_tss_source(tss_source, genome)
  chr_len_db <- get("CHR_LEN_DB", envir = asNamespace("rGREAT"))
  if (!src$genome %in% names(chr_len_db)) stop("rGREAT has no chromosome lengths for ", src$genome, call. = FALSE)
  chrlen <- chr_len_db[[src$genome]]

  if (src$category == "GREAT") {
    f <- system.file("extdata", paste0("GREATv4.genes.", src$genome, ".tsv.gz"), package = "rGREAT")
    if (!nzchar(f)) stop("GREAT v4 gene table for ", src$genome, " not found in rGREAT.", call. = FALSE)
    tb <- utils::read.table(gzfile(f), stringsAsFactors = FALSE, quote = "", comment.char = "#")
    tb[is.na(tb[, 1]), 1] <- ""
    entrez <- .great_map_table_to_entrez(tb, src$genome)
    gr <- GenomicRanges::GRanges(tb[, 2], IRanges::IRanges(tb[, 3] + 1, tb[, 3] + 1),
                                 strand = tb[, 4], gene_id = entrez, symbol = tb[, 5])
    n_all <- length(gr)
    gr <- gr[!is.na(gr$gene_id)]
    gr <- gr[as.character(GenomicRanges::seqnames(gr)) %in% names(chrlen)]
    gr <- gr[!duplicated(gr$gene_id)]
    if (verbose) message("[GREAT] GREAT v4 ", src$genome, " table: ", n_all, " genes; ",
                         length(gr), " with Entrez IDs on standard chromosomes")
    GenomeInfoDb::seqlevels(gr) <- names(chrlen)
    GenomeInfoDb::seqlengths(gr) <- chrlen
    symbol_map <- stats::setNames(gr$symbol, gr$gene_id)
    S4Vectors::mcols(gr) <- S4Vectors::mcols(gr)[, "gene_id", drop = FALSE]
    ext <- rGREAT::extendTSS(gr, seqlengths = chrlen, genome = src$genome, gene_id_type = "ENTREZ",
                             mode = mode, basal_upstream = basal_upstream,
                             basal_downstream = basal_downstream, extension = extension,
                             verbose = FALSE)
  } else {
    if (!requireNamespace(src$txdb_pkg, quietly = TRUE)) stop(src$txdb_pkg, " is not installed.", call. = FALSE)
    ext <- rGREAT::extendTSSFromTxDb(src$txdb_pkg, verbose = FALSE, mode = mode,
                                     basal_upstream = basal_upstream,
                                     basal_downstream = basal_downstream, extension = extension)
    info <- gene_info_from_orgdb(unique(ext$gene_id), src$genome)
    symbol_map <- stats::setNames(info$symbol_orgdb, info$gene_entrezID)
    if (verbose) message("[GREAT] ", src$txdb_pkg, ": ", length(unique(ext$gene_id)), " protein-coding genes")
  }
  ext$gene_id <- as.character(ext$gene_id)
  attr(ext, "symbol_map") <- symbol_map
  attr(ext, "great_params") <- list(tss_source = src$label, genome = src$genome, mode = mode,
                                    basal_upstream = basal_upstream, basal_downstream = basal_downstream,
                                    extension = extension,
                                    rGREAT_version = as.character(utils::packageVersion("rGREAT")),
                                    n_genes = length(unique(ext$gene_id)))
  ext
}

# ------------------------------------------------------------
# Region -> gene associations (long format), same rule as
# rGREAT::getRegionGeneAssociations(): full-region overlap with each gene's
# regulatory domain; signed distance from the nearest region edge to TSS.
# ------------------------------------------------------------
great_region_gene_associations <- function(regions_df, extended_tss) {
  gr <- GenomicRanges::GRanges(regions_df$chr, IRanges::IRanges(regions_df$start, regions_df$end))
  ov <- GenomicRanges::findOverlaps(gr, extended_tss, ignore.strand = TRUE)
  qh <- S4Vectors::queryHits(ov); sh <- S4Vectors::subjectHits(ov)
  st <- regions_df$start[qh]; en <- regions_df$end[qh]
  tp <- extended_tss$tss_position[sh]; ts <- as.character(extended_tss$tss_strand[sh])
  up   <- (en < tp & ts == "+") | (st > tp & ts == "-")
  down <- (en < tp & ts == "-") | (st > tp & ts == "+")
  dist <- pmin(abs(st - tp), abs(en - tp))
  dist[up] <- -dist[up]; dist[!(up | down)] <- 0L
  sym <- attr(extended_tss, "symbol_map")
  gid <- extended_tss$gene_id[sh]
  out <- data.frame(RegionID = regions_df$RegionID[qh], chr = regions_df$chr[qh],
                    start = st, end = en, gene_entrezID = gid,
                    gene_symbol = unname(ifelse(is.na(sym[gid]), gid, sym[gid])),
                    distance_to_TSS = as.integer(dist), stringsAsFactors = FALSE)
  if ("module" %in% colnames(regions_df)) out$module <- regions_df$module[qh]
  out[order(match(out$RegionID, regions_df$RegionID), abs(out$distance_to_TSS)), , drop = FALSE]
}

# Collapse long associations to one row per region, comethyl style
# ("A | B"), ordered by absolute distance to TSS (closest first).
great_collapse_associations <- function(assoc, regions_df, gene_info = NULL) {
  key <- assoc$RegionID
  col <- function(v) tapply(v, factor(key, levels = unique(key)), function(x) paste(x, collapse = " | "))
  agg <- data.frame(RegionID = unique(key),
                    gene_symbol = as.character(col(assoc$gene_symbol)),
                    distance_to_TSS = as.character(col(assoc$distance_to_TSS)),
                    gene_entrezID = as.character(col(assoc$gene_entrezID)),
                    n_genes = as.integer(tapply(assoc$gene_symbol, factor(key, levels = unique(key)), length)),
                    stringsAsFactors = FALSE)
  if (!is.null(gene_info) && nrow(gene_info) > 0) {
    gi <- gene_info[match(assoc$gene_entrezID, gene_info$gene_entrezID), ]
    agg$gene_description <- as.character(col(ifelse(is.na(gi$gene_description), "", gi$gene_description)))
    agg$gene_ensemblID   <- as.character(col(ifelse(is.na(gi$gene_ensemblID), "", gi$gene_ensemblID)))
  } else {
    agg$gene_description <- NA_character_; agg$gene_ensemblID <- NA_character_
  }
  out <- data.frame(RegionID = regions_df$RegionID, stringsAsFactors = FALSE)
  out <- merge(out, agg, by = "RegionID", all.x = TRUE, sort = FALSE)
  out <- out[match(regions_df$RegionID, out$RegionID), , drop = FALSE]
  out$n_genes[is.na(out$n_genes)] <- 0L
  rownames(out) <- NULL
  out
}

# ================================================================
# GENE-SET CACHE (GO BP/MF/CC, KEGG, Reactome, MSigDB Hallmark)
#
# One RDS holding every collection as Entrez-ID gene sets, plus term
# names, a symbol<->Entrez map and provenance. Built once on a node with
# internet (KEGG and MSigDB downloads), then loaded by every 12b run so
# GREAT, clusterProfiler and the gene universe all use identical,
# versioned gene sets.
#
# source = "native": GO from org.Hs.eg.db (GOALL; same data
#            clusterProfiler::enrichGO uses), KEGG from the KEGG REST API
#            (clusterProfiler::download_KEGG), Reactome from reactome.db
#            (as ReactomePA), Hallmark from msigdbr.
# source = "msigdb": everything from one MSigDB release via msigdbr
#            (C5 GO:BP/MF/CC, C2 CP:KEGG_MEDICUS, C2 CP:REACTOME, H).
# ================================================================
GENE_SET_COLLECTIONS <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")

# ------------------------------------------------------------
# parse_collections
#   Validate a comma list of collections (case-insensitive) against
#   GENE_SET_COLLECTIONS; NULL, "" or "all" -> all six.
# ------------------------------------------------------------
parse_collections <- function(x) {
  if (is.null(x) || !nzchar(x) || tolower(x) == "all") return(GENE_SET_COLLECTIONS)
  v <- trimws(strsplit(x, ",")[[1]])
  m <- GENE_SET_COLLECTIONS[match(tolower(v), tolower(GENE_SET_COLLECTIONS))]
  if (anyNA(m)) stop("Unknown collection(s): ", paste(v[is.na(m)], collapse = ", "),
                     ". Choose from: ", paste(GENE_SET_COLLECTIONS, collapse = ", "), call. = FALSE)
  unique(m)
}

# msigdbr argument names changed in v10 (category -> collection,
# subcategory -> subcollection); support both.
.msigdbr_fetch <- function(collection, subcollection = NULL) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr is not installed.", call. = FALSE)
  fa <- names(formals(msigdbr::msigdbr))
  args <- list(species = "Homo sapiens")
  if ("collection" %in% fa) { args$collection <- collection; if (!is.null(subcollection)) args$subcollection <- subcollection }
  else { args$category <- collection; if (!is.null(subcollection)) args$subcategory <- subcollection }
  d <- do.call(msigdbr::msigdbr, args)
  idcol <- intersect(c("ncbi_gene", "entrez_gene", "human_entrez_gene"), colnames(d))[1]
  symcol <- intersect(c("gene_symbol", "human_gene_symbol"), colnames(d))[1]
  data.frame(term_id = d$gs_name, term_name = d$gs_name, entrez = as.character(d[[idcol]]),
             symbol = d[[symcol]], stringsAsFactors = FALSE)
}

# ------------------------------------------------------------
# .long_to_sets
#   Internal: long table (term_id, term_name, entrez[, symbol]) ->
#   list(sets = term -> unique Entrez IDs, names = term -> name,
#   symbol_pairs = unique symbol/Entrez pairs if present).
# ------------------------------------------------------------
.long_to_sets <- function(df) {
  df <- df[!is.na(df$entrez) & nzchar(df$entrez), ]
  sets <- lapply(split(df$entrez, df$term_id), unique)
  nm <- df$term_name[!duplicated(df$term_id)]; names(nm) <- df$term_id[!duplicated(df$term_id)]
  out <- list(sets = sets, names = nm[names(sets)])
  if ("symbol" %in% colnames(df)) out$symbol_pairs <- unique(df[!is.na(df$symbol), c("symbol", "entrez")])
  out
}

# ------------------------------------------------------------
# .native_collection
#   Internal: build one collection from its native source -- GO
#   BP/MF/CC from org.Hs.eg.db GOALL (+ GO.db names), KEGG from the KEGG
#   REST API (internet), Reactome from reactome.db (R-HSA only), Hallmark
#   from msigdbr. Returns .long_to_sets() output.
# ------------------------------------------------------------
.native_collection <- function(coll) {
  if (coll %in% c("GO_BP", "GO_MF", "GO_CC")) {
    for (p in c("org.Hs.eg.db", "GO.db")) if (!requireNamespace(p, quietly = TRUE)) stop(p, " is required for GO.", call. = FALSE)
    ont <- sub("GO_", "", coll)
    m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
           keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "GOALL"),
           keytype = "GOALL", columns = c("ENTREZID", "ONTOLOGYALL")))
    m <- m[m$ONTOLOGYALL == ont & !is.na(m$ENTREZID), ]
    terms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(m$GOALL),
                                                    keytype = "GOID", columns = "TERM"))
    df <- data.frame(term_id = m$GOALL, entrez = m$ENTREZID, stringsAsFactors = FALSE)
    df$term_name <- terms$TERM[match(df$term_id, terms$GOID)]
    return(.long_to_sets(df))
  }
  if (coll == "KEGG") {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("clusterProfiler is required for KEGG.", call. = FALSE)
    kg <- clusterProfiler::download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")
    df <- data.frame(term_id = kg$KEGGPATHID2EXTID[[1]], entrez = as.character(kg$KEGGPATHID2EXTID[[2]]),
                     stringsAsFactors = FALSE)
    df$term_name <- kg$KEGGPATHID2NAME[[2]][match(df$term_id, kg$KEGGPATHID2NAME[[1]])]
    return(.long_to_sets(df))
  }
  if (coll == "Reactome") {
    if (!requireNamespace("reactome.db", quietly = TRUE)) stop("reactome.db is required for Reactome.", call. = FALSE)
    p2g <- as.list(reactome.db::reactomePATHID2EXTID)
    p2g <- p2g[grepl("^R-HSA-", names(p2g))]
    nm <- unlist(as.list(reactome.db::reactomePATHID2NAME)[names(p2g)])
    nm <- sub("^Homo sapiens: ", "", nm)
    sets <- lapply(p2g, function(x) unique(as.character(x)))
    return(list(sets = sets, names = nm[names(sets)]))
  }
  if (coll == "Hallmark") return(.long_to_sets(.msigdbr_fetch("H")))
  stop("Unknown collection: ", coll, call. = FALSE)
}

# ------------------------------------------------------------
# .msigdb_collection
#   Internal: build one collection from MSigDB via msigdbr (C5 GO,
#   C2 CP:KEGG_MEDICUS, C2 CP:REACTOME, H).
# ------------------------------------------------------------
.msigdb_collection <- function(coll) {
  spec <- switch(coll,
    GO_BP = list("C5", "GO:BP"), GO_MF = list("C5", "GO:MF"), GO_CC = list("C5", "GO:CC"),
    KEGG = list("C2", "CP:KEGG_MEDICUS"), Reactome = list("C2", "CP:REACTOME"), Hallmark = list("H", NULL))
  .long_to_sets(.msigdbr_fetch(spec[[1]], spec[[2]]))
}

# ------------------------------------------------------------
# build_gene_set_cache
#   Build every requested collection (source = "native" or
#   "msigdb"), add a genome-wide symbol -> Entrez map (org.Hs.eg.db, else
#   MSigDB pairs), package versions and build time, and saveRDS() it to
#   cache_file. Run once on a node with internet. Returns it invisibly.
# ------------------------------------------------------------
build_gene_set_cache <- function(collections = GENE_SET_COLLECTIONS, source = c("native", "msigdb"),
                                 cache_file, verbose = TRUE) {
  source <- match.arg(source)
  out <- list(collections = list(), term_names = list(), meta = list(
    source = source, built = timestamp_now(), collections = collections))
  pairs <- NULL
  for (coll in collections) {
    if (verbose) message("[gene_sets] Building ", coll, " (", source, ")")
    x <- if (source == "native") .native_collection(coll) else .msigdb_collection(coll)
    out$collections[[coll]] <- x$sets; out$term_names[[coll]] <- x$names
    if (!is.null(x$symbol_pairs)) pairs <- unique(rbind(pairs, x$symbol_pairs))
    if (verbose) message("[gene_sets]   ", length(x$sets), " sets")
  }
  vers <- c("org.Hs.eg.db", "GO.db", "reactome.db", "msigdbr", "clusterProfiler", "KEGGREST")
  out$meta$package_versions <- vapply(vers, function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, "")
  # symbol -> Entrez map for converting gene lists (symbols) to Entrez:
  # org.Hs.eg.db (genome-wide) if available, else the pairs seen in MSigDB.
  sm <- tryCatch({
    m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
           keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "ENTREZID"),
           keytype = "ENTREZID", columns = "SYMBOL"))
    out$meta$symbol_map_source <- "org.Hs.eg.db"
    stats::setNames(m$ENTREZID, m$SYMBOL)
  }, error = function(e) NULL)
  if (is.null(sm) && !is.null(pairs)) {
    pairs <- pairs[!duplicated(pairs$symbol), ]
    sm <- stats::setNames(pairs$entrez, pairs$symbol)
    out$meta$symbol_map_source <- "MSigDB gene tables (org.Hs.eg.db unavailable)"
    if (verbose) message("[gene_sets] WARNING: org.Hs.eg.db unavailable; symbol map built from MSigDB gene tables")
  }
  out$symbol_map <- sm
  safe_dir_create(dirname(cache_file))
  saveRDS(out, cache_file)
  if (verbose) message("[gene_sets] Saved: ", cache_file)
  invisible(out)
}

# ------------------------------------------------------------
# load_gene_set_cache
#   readRDS() a cache from build_gene_set_cache(); stops if it isn't
#   one or lacks any requested collection.
# ------------------------------------------------------------
load_gene_set_cache <- function(cache_file, collections) {
  validate_file_exists(cache_file, "--gene_set_cache")
  gs <- readRDS(cache_file)
  if (!all(c("collections", "term_names", "meta") %in% names(gs))) {
    stop("Not a gene-set cache built by build_gene_set_cache(): ", cache_file, call. = FALSE)
  }
  miss <- setdiff(collections, names(gs$collections))
  if (length(miss) > 0) stop("Gene-set cache lacks: ", paste(miss, collapse = ", "),
                             ". Rebuild it with those collections.", call. = FALSE)
  gs
}

# ------------------------------------------------------------
# symbols_to_entrez
#   Map gene symbols to Entrez via the cache's symbol map. Returns
#   list(entrez = unique mapped IDs, n_in, n_mapped). Symbols with no
#   mapping (e.g. outdated names) are dropped silently -- check n_mapped.
# ------------------------------------------------------------
symbols_to_entrez <- function(symbols, gene_set_cache) {
  symbols <- unique(trimws(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (is.null(gene_set_cache$symbol_map)) stop("Gene-set cache has no symbol map (org.Hs.eg.db was missing at build time).", call. = FALSE)
  e <- unname(gene_set_cache$symbol_map[symbols])
  list(entrez = unique(e[!is.na(e)]), n_in = length(symbols), n_mapped = sum(!is.na(e)))
}

# ================================================================
# FAST GREAT TEST (same statistics as rGREAT::great(), all modules at once)
#
# rGREAT::great() recomputes every gene set's domain coverage of the
# background on each call (~7 min per module for ~13k gene sets and a
# 465k-region background). That coverage does not depend on the module,
# so great_fast_prepare() computes it once and great_fast_test() then
# only counts hits per module with sparse matrices. The procedure is
# rGREAT 2.4's, step by step:
#   background  = reduce(background); total = sum(width)
#   domains     = extended_tss clipped to background (pintersect)
#   gene sets   = restricted to genes with clipped domains, then
#                 min_gene_set_size filter (gene_set_size reported after this)
#   prop        = width(reduce(clipped domains of set genes)) / total
#   regions     = midpoints, intersect()-ed with background (n_total)
#   n_hits      = regions overlapping the set's clipped domains
#   p           = 1 - pbinom(n_hits - 1, n_total, prop); fold = n_hits / (n_total * prop)
#   gene hits   = set genes whose domains contain >= 1 region
#   p_hyper     = 1 - phyper(gene_hits - 1, set_size, n_gene_total - set_size, n_gene_gr)
# 12b's --great_engine rgreat runs rGREAT::great() itself for verification.
# ================================================================
great_fast_prepare <- function(extended_tss, background_gr, gene_sets, min_size = 5, max_size = Inf,
                               verbose = TRUE) {
  bg <- GenomicRanges::reduce(GenomicRanges::granges(background_gr), ignore.strand = TRUE)
  GenomicRanges::strand(bg) <- "*"
  total <- sum(as.numeric(GenomicRanges::width(bg)))
  ext <- extended_tss; GenomicRanges::strand(ext) <- "*"
  ov <- GenomicRanges::findOverlaps(ext, bg)
  pieces <- GenomicRanges::pintersect(ext[S4Vectors::queryHits(ov)], bg[S4Vectors::subjectHits(ov)])
  pieces$gene_id <- ext$gene_id[S4Vectors::queryHits(ov)]
  pieces$tss_position <- ext$tss_position[S4Vectors::queryHits(ov)]
  genes <- sort(unique(pieces$gene_id))

  sets <- lapply(gene_sets, function(x) intersect(x, genes))
  n <- lengths(sets)
  sets <- sets[n >= min_size & n <= max_size]
  if (length(sets) == 0) stop("No gene set left after size filtering.", call. = FALSE)
  G <- Matrix::sparseMatrix(i = match(unlist(sets, use.names = FALSE), genes),
                            j = rep(seq_along(sets), lengths(sets)), x = 1,
                            dims = c(length(genes), length(sets)), dimnames = list(genes, names(sets)))

  # Disjoint segments of all clipped domains, grouped by the exact set of
  # genes covering them ("combo"); most combos are a single gene.
  D <- GenomicRanges::disjoin(pieces, ignore.strand = TRUE)
  dov <- GenomicRanges::findOverlaps(D, pieces)
  seg_genes <- split(pieces$gene_id[S4Vectors::subjectHits(dov)], S4Vectors::queryHits(dov))
  key <- vapply(seg_genes, function(g) paste(sort(unique(g)), collapse = ";"), "")
  seg_key <- character(length(D)); seg_key[as.integer(names(seg_genes))] <- key
  combos <- unique(key)
  seg_combo <- match(seg_key, combos)
  combo_width <- as.numeric(tapply(as.numeric(GenomicRanges::width(D)), factor(seg_combo, levels = seq_along(combos)), sum))
  combo_genes <- strsplit(combos, ";", fixed = TRUE)
  K <- Matrix::sparseMatrix(i = rep(seq_along(combos), lengths(combo_genes)),
                            j = match(unlist(combo_genes), genes), x = 1,
                            dims = c(length(combos), length(genes)))
  Pbin <- K %*% G; Pbin@x[] <- 1                      # combo covered by set (any gene)
  prop <- as.numeric(Matrix::crossprod(Pbin, combo_width)) / total
  if (verbose) message("[GREAT-fast] background ", length(bg), " ranges (", round(total / 1e6, 2), " Mb); ",
                       length(genes), " genes with domains in background; ", length(sets), " gene sets; ",
                       length(combos), " domain combos")
  list(bg = bg, total = total, pieces = pieces, genes = genes, sets = sets, G = G, D = D,
       seg_combo = seg_combo, K = K, Pbin = Pbin, prop = prop)
}

# ------------------------------------------------------------
# great_fast_test
#   Run the GREAT tests for every module at once using a
#   great_fast_prepare() object. regions_by_module = named list of region
#   data.frames (chr, start, end). Region midpoints are used, as in rGREAT.
#   Returns a named list of per-module tables (binomial region test +
#   hypergeometric gene test, BH within module) with attributes
#   "hit_genes" (Entrez IDs hit) and "n_total" (midpoints in background).
# ------------------------------------------------------------
great_fast_test <- function(prep, regions_by_module, verbose = TRUE) {
  mods <- names(regions_by_module)
  nm <- length(mods); nc <- nrow(prep$K); ng <- length(prep$genes)
  Cm <- matrix(0, nc, nm); Sg <- matrix(0, ng, nm); Ng <- matrix(0, ng, nm)
  n_total <- n_gene_gr <- numeric(nm)
  for (j in seq_len(nm)) {
    r <- regions_by_module[[j]]
    mid <- GenomicRanges::GRanges(r$chr, IRanges::IRanges(floor((r$start + r$end) / 2), width = 1))
    m <- GenomicRanges::intersect(mid, prep$bg, ignore.strand = TRUE)       # as rGREAT
    n_total[j] <- length(m)
    so <- GenomicRanges::findOverlaps(m, prep$D, select = "first")
    cc <- prep$seg_combo[so[!is.na(so)]]
    tb <- tabulate(cc, nbins = nc); Cm[, j] <- tb
    po <- GenomicRanges::findOverlaps(m, prep$pieces)
    gi <- match(prep$pieces$gene_id[S4Vectors::subjectHits(po)], prep$genes)
    d <- abs(prep$pieces$tss_position[S4Vectors::subjectHits(po)] - GenomicRanges::start(m)[S4Vectors::queryHits(po)])
    Sg[, j] <- as.numeric(tapply(d, factor(gi, levels = seq_len(ng)), sum, default = 0))
    Ng[, j] <- tabulate(gi, ng)
    n_gene_gr[j] <- sum(Ng[, j] > 0)
  }
  hits <- as.matrix(Matrix::crossprod(prep$Pbin, Matrix::Matrix(Cm, sparse = TRUE)))     # sets x modules
  gene_hit <- Matrix::Matrix((Ng > 0) * 1, sparse = TRUE)
  ghits <- as.matrix(Matrix::crossprod(prep$G, gene_hit))
  sumd <- as.matrix(Matrix::crossprod(prep$G, Matrix::Matrix(Sg, sparse = TRUE)))
  nd <- as.matrix(Matrix::crossprod(prep$G, Matrix::Matrix(Ng, sparse = TRUE)))
  size <- lengths(prep$sets); n_gene_total <- length(prep$genes)
  out <- list()
  for (j in seq_len(nm)) {
    n_hits <- hits[, j]
    p <- ifelse(n_hits == 0, 1, 1 - stats::pbinom(n_hits - 1, n_total[j], prep$prop))
    df <- data.frame(id = names(prep$sets), genome_fraction = prep$prop, observed_region_hits = n_hits,
                     fold_enrichment = n_hits / (prep$prop * n_total[j]), p_value = p,
                     p_adjust = stats::p.adjust(p, "BH"),
                     mean_tss_dist = round(ifelse(nd[, j] > 0, sumd[, j] / nd[, j], NaN)),
                     observed_gene_hits = ghits[, j], gene_set_size = size, stringsAsFactors = FALSE)
    df$fold_enrichment_hyper <- df$observed_gene_hits / (df$gene_set_size * n_gene_gr[j] / n_gene_total)
    df$p_value_hyper <- 1 - stats::phyper(df$observed_gene_hits - 1, df$gene_set_size,
                                          n_gene_total - df$gene_set_size, n_gene_gr[j])
    df$p_adjust_hyper <- stats::p.adjust(df$p_value_hyper, "BH")
    rownames(df) <- NULL
    attr(df, "hit_genes") <- prep$genes[Ng[, j] > 0]
    attr(df, "n_total") <- n_total[j]
    out[[mods[j]]] <- df
  }
  out
}