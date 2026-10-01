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
#
# NOTES
#   - These helpers are generic and do not include project-specific
#     trait recoding rules.
#   - sample_info is assumed to be analysis-ready before entering the
#     pipeline.
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

  # First Excel column contains your sample IDs / original rownames
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
# Returns a list with:
#   requested
#   found
#   missing
# ------------------------------------------------------------
#old before trait files with dots and spaces
# resolveTraits <- function(requested_traits,
#                           available_traits,
#                           label = "traits",
#                           verbose = TRUE) {
#   if (is.null(requested_traits)) {
#     return(list(
#       requested = character(0),
#       found = character(0),
#       missing = character(0)
#     ))
#   }

#   requested_traits <- unique(trimws(requested_traits))
#   requested_traits <- requested_traits[requested_traits != ""]

#   found <- intersect(requested_traits, available_traits)
#   missing <- setdiff(requested_traits, available_traits)

#   if (verbose) {
#     message("[resolveTraits] ", label, ": requested = ", length(requested_traits),
#             ", found = ", length(found),
#             ", missing = ", length(missing))
#   }

#   return(list(
#     requested = requested_traits,
#     found = found,
#     missing = missing
#   ))
# }

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

  # Create a comparison key that ignores differences such as:
  # spaces vs dots, parentheses, colon, slash, hyphen, etc.
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
                           axis.text.size = 18, legend.position = c(1.08, 0.915),
                           legend.text.size = 18, legend.title.size = 18,
                           colColorMargins = c(-0.7,4.21,1.2,11.07),
                           save = TRUE,
                           file = "ME_Trait_Correlation_Heatmap.pdf",
                           width = 11, height = 9.5, verbose = TRUE,
                           ## NEW controls for the labeled color bar
                           showColorBar = TRUE,
                           showColorBarLabels = TRUE,             # add names on boxes
                           colorBarLabelPos = c("inside","below"),
                           colorBarLabelSize = 4,                 # font size of the names on the color boxes
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
  # Hide axis labels if we’re labeling the color bar (to avoid duplicates) #i changed this if (showColorBar && showColorBarLabels) to if (showColorBar) to show the color bar even if the labels are not shown
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

# #below are some of the settings that can be used to make plots
# plotMEtraitCor(
#   MEtraitCor,                                          # data frame with columns: module, trait, p, and cor/bicor

#   moduleOrder       = 1:length(unique(MEtraitCor$module)), # order of module columns (left→right)
#   traitOrder        = 1:length(unique(MEtraitCor$trait)),  # order of trait rows (bottom→top; reversed internally)

#   topOnly           = FALSE,                           # if TRUE, plot only the nTop most significant cells
#   nTop              = 15,                              # number of most-significant cells when topOnly=TRUE
#   p                 = 0.05,                            # significance cutoff used for stars/p-values overlay

#   label.type        = c("star", "p"),                  # overlay type; default resolves to "star"
#   label.size        = 8,                               # size of the star or p-value text
#   label.nudge_y     = -0.38,                           # vertical nudge for label positioning

#   colors            = blueWhiteRed(100, gamma = 0.9),  # heatmap palette (min→white→max)
#   limit             = NULL,                            # color scale limit; NULL = max(abs(cor)) auto

#   axis.text.size    = 12,                              # font size of axis tick labels (modules/traits)
#   legend.position   = c(1.08, 0.915),                  # legend position (inside plotting area)
#   legend.text.size  = 12,                              # font size of legend tick labels
#   legend.title.size = 16,                              # font size of legend title

#   colColorMargins   = c(-0.7, 4.21, 1.2, 11.07),       # margins (lines) around the module color bar (t,r,b,l)
#   save              = TRUE,                            # write the figure to disk
#   file              = "ME_Trait_Correlation_Heatmap.pdf", # output filename (when save=TRUE)
#   width             = 11,                              # output width in inches
#   height            = 9.5,                             # output height in inches
#   verbose           = TRUE,                            # print progress messages

#   # ---- NEW labeled color-bar controls ----
#   showColorBar        = TRUE,                          # draw the module color strip under the heatmap
#   showColorBarLabels  = TRUE,                          # put module names on the color boxes
#   colorBarLabelPos    = c("inside","below"),           # where to place names; default resolves to "inside"
#   colorBarLabelSize   = 3.2,                           # font size of the names on the color boxes
#   colorBarLabelAngle  = 90,                            # rotation of those names (90 = vertical)
#   colorBarRelHeight   = 0.10,                          # relative height of the color bar vs heatmap
#   syncWidths          = TRUE,                          # force exact column alignment via gtable widths
#   autoContrastLabels  = TRUE                           # auto-pick black/white text for readability on each color
# )


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


# #below are some of the settings that can be used to make plots
# plotMEtraitCor(
#   MEtraitCor,                                          # data frame with columns: module, trait, p, and cor/bicor

#   moduleOrder       = 1:length(unique(MEtraitCor$module)), # order of module columns (left→right)
#   traitOrder        = 1:length(unique(MEtraitCor$trait)),  # order of trait rows (bottom→top; reversed internally)

#   topOnly           = FALSE,                           # if TRUE, plot only the nTop most significant cells
#   nTop              = 15,                              # number of most-significant cells when topOnly=TRUE
#   p                 = 0.05,                            # significance cutoff used for stars/p-values overlay

#   label.type        = c("star", "p"),                  # overlay type; default resolves to "star"
#   label.size        = 8,                               # size of the star or p-value text
#   label.nudge_y     = -0.38,                           # vertical nudge for label positioning

#   colors            = blueWhiteRed(100, gamma = 0.9),  # heatmap palette (min→white→max)
#   limit             = NULL,                            # color scale limit; NULL = max(abs(cor)) auto

#   axis.text.size    = 12,                              # font size of axis tick labels (modules/traits)
#   legend.position   = c(1.08, 0.915),                  # legend position (inside plotting area)
#   legend.text.size  = 12,                              # font size of legend tick labels
#   legend.title.size = 16,                              # font size of legend title

#   colColorMargins   = c(-0.7, 4.21, 1.2, 11.07),       # margins (lines) around the module color bar (t,r,b,l)
#   save              = TRUE,                            # write the figure to disk
#   file              = "ME_Trait_Correlation_Heatmap.pdf", # output filename (when save=TRUE)
#   width             = 11,                              # output width in inches
#   height            = 9.5,                             # output height in inches
#   verbose           = TRUE,                            # print progress messages

#   # ---- NEW labeled color-bar controls ----
#   showColorBar        = TRUE,                          # draw the module color strip under the heatmap
#   showColorBarLabels  = TRUE,                          # put module names on the color boxes
#   colorBarLabelPos    = c("inside","below"),           # where to place names; default resolves to "inside"
#   colorBarLabelSize   = 3.2,                           # font size of the names on the color boxes
#   colorBarLabelAngle  = 90,                            # rotation of those names (90 = vertical)
#   colorBarRelHeight   = 0.10,                          # relative height of the color bar vs heatmap
#   syncWidths          = TRUE,                          # force exact column alignment via gtable widths
#   autoContrastLabels  = TRUE                           # auto-pick black/white text for readability on each color
# )


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