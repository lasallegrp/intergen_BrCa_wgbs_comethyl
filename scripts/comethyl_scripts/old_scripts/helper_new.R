# ================================================================
# helper.R
#
# Shared helper functions for the reproducible comethyl pipeline.
#
# FUNCTIONS INCLUDED
#   - readTraitFile()
#   - readSampleInfo()
#   - resolveTraits()
#   - plotPCTrait()
#   - plotTraitDendrogramFromPC()
#   - plotMEtraitCor()
#   - save_me_trait_method_outputs()
#   - plotMEtraitViolin() and related ME-trait pair helpers
#   - collect_set_files() / read_trait_set_file() / get_set_name()
#   - ANNOTATION HELPERS (region -> gene symbol):
#       get_arg(), trim_or_null(), safe_dir_create(), timestamp_now(),
#       append_log(), write_lines_safe(), stop_if_missing(),
#       validate_file_exists(), validate_regions_df(),
#       offline_nearest_gene(), annotate_offline_only(),
#       annotate_regions_safe()
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
# ================================================================

suppressPackageStartupMessages({
  library(openxlsx)
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
    df <- openxlsx::read.xlsx(file, rowNames = TRUE)
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

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

  requested_traits <- unique(trimws(requested_traits))
  requested_traits <- requested_traits[requested_traits != ""]

  found <- intersect(requested_traits, available_traits)
  missing <- setdiff(requested_traits, available_traits)

  if (verbose) {
    message("[resolveTraits] ", label, ": requested = ", length(requested_traits),
            ", found = ", length(found),
            ", missing = ", length(missing))
  }

  return(list(
    requested = requested_traits,
    found = found,
    missing = missing
  ))
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
# Generalized heatmap for PC-trait association statistics.
#
# MODES
#   1) selected_traits
#      - plot user-supplied traits only
#
#   2) top_cells
#      - keep top_n_cells PC-trait pairs by smallest p-value
#
#   3) top_traits
#      - rank traits by smallest p-value across all PCs
#      - keep top_n_traits traits
#      - plot all PCs against those traits
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
    stop("[plotMEtraitCor_pub] Need 'cor' or 'bicor' column")

  if (verbose) message("[plotMEtraitCor_pub] corType = ", corType)

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
    if (verbose) message(sprintf("[plotMEtraitCor_pub] autoSize -> width=%.2f, height=%.2f (nX=%d,nY=%d)",
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

parse_bool <- function(x, arg_name = "argument") {
  if (is.logical(x)) return(x)
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE. Got: ", x)
}

write_log_lines <- function(lines, file) {
  writeLines(as.character(lines), con = file)
}

write_vector_file <- function(x, file) {
  x <- unique(as.character(x))
  if (length(x) == 0) x <- character(0)
  writeLines(x, con = file)
}

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
# Builds a ME-trait heatmap and can adjust colors 
# ------------------------------------------------------------
plotMEtraitCor <- function(MEtraitCor,
                           moduleOrder = 1:length(unique(MEtraitCor$module)),
                           traitOrder  = 1:length(unique(MEtraitCor$trait)),
                           topOnly = FALSE, nTop = 15, p = 0.05,
                           label.type = c("star", "p"), label.size = 8,
                           label.nudge_y = -0.38,
                           colors = blueWhiteRed(100, gamma = 0.9), limit = NULL,
                           axis.text.size = 12, legend.position = c(1.08, 0.915),
                           legend.text.size = 12, legend.title.size = 16,
                           colColorMargins = c(-0.7,4.21,1.2,11.07),
                           save = TRUE,
                           file = "ME_Trait_Correlation_Heatmap.pdf",
                           width = 11, height = 9.5, verbose = TRUE,
                           ## NEW controls for the labeled color bar
                           showColorBar = TRUE,
                           showColorBarLabels = TRUE,             # add names on boxes
                           colorBarLabelPos = c("inside","below"),
                           colorBarLabelSize = 3.2,
                           colorBarLabelAngle = 90,                # 90 = vertical
                           colorBarRelHeight = 0.10,               # bar height rel to heatmap
                           syncWidths = TRUE,                      # force exact alignment
                           autoContrastLabels = TRUE               # black/white text auto
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
  # Hide axis labels if we’re labeling the color bar (to avoid duplicates)
  if (showColorBar && showColorBarLabels) {
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

auto_ylim <- function(ME_vec) {
  y <- ME_vec[is.finite(ME_vec)]
  if (length(y) == 0) return(c(-1, 1))
  r <- range(y)
  pad <- 0.08 * diff(r)
  if (!is.finite(pad) || pad == 0) pad <- 0.1
  c(r[1] - pad, r[2] + pad)
}

is_binary_like <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) == 2
}

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

read_trait_set_file <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

get_set_name <- function(file) {
  tools::file_path_sans_ext(basename(file))
}

# ================================================================
# ANNOTATION HELPERS (region -> gene symbol, genomic location, CpG context)
#
# Added for scripts that annotate comethyl regions with nearest
# gene symbols, e.g. 12a_annotate_modules.R (module-level regions)
# and 12c_annotate_region_matrix.R (all filtered regions,
# independent of module assignment).
#
# annotate_offline_only() / annotate_regions_safe() now also add,
# when possible:
#   genomic_location  Promoter / 5UTR / 3UTR / Exon / Intron /
#                     Downstream / Intergenic, via ChIPseeker +
#                     a local TxDb (no internet needed at call time,
#                     only for the one-time package install).
#   cpg_context       Island / Shore / Shelf / OpenSea, via a local
#                     UCSC CpG island table. Download once from ANY
#                     internet-connected machine (your laptop is
#                     fine) and pass its path as cpg_island_file:
#                       https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz
#                     If cpg_island_file is NULL, this column is
#                     skipped (no error).
#
# Calls GenomicRanges/IRanges/GenomeInfoDb/S4Vectors/AnnotationDbi/
# ChIPseeker/comethyl functions via pkg::fn() namespacing throughout,
# so those packages only need to be installed, not attached, in
# whatever script sources this file.
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

file_checksum <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

# ------------------------------------------------------------
# CLI argument helpers
# ------------------------------------------------------------
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

trim_or_null <- function(x) {
  if (is.null(x) || is.na(x)) return(NULL)
  x <- trimws(x)
  if (!nzchar(x)) return(NULL)
  x
}

# ------------------------------------------------------------
# Filesystem / logging helpers
# ------------------------------------------------------------
safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

append_log <- function(logfile, ...) {
  txt <- paste0("[", timestamp_now(), "] ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = logfile, append = TRUE)
}

write_lines_safe <- function(x, file) {
  writeLines(as.character(x), con = file, useBytes = TRUE)
}

# ------------------------------------------------------------
# Validation helpers
# ------------------------------------------------------------
stop_if_missing <- function(x, label) {
  if (is.null(x) || !nzchar(x)) stop("Missing required argument: ", label, call. = FALSE)
}

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
# Annotation core
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