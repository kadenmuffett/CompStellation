#' Plot a Star Plot of Group Ordination Position Surrounded by a Taxonomic Iris Plot
#'
#' This experimental (beta) function combines an ordination radar plot (the "star")
#' in the center, and surrounds it with a circular stacked bar chart (the "iris plot")
#' representing the taxonomic composition of all individual samples belonging to that group.
#'
#' @section Interpreting the plot:
#' The two regions of the radial axis carry different units: the centre shows
#' ordination scores (offset so that negative scores are plottable, with the
#' tick labels undoing the offset) and the outer ring shows relative abundance
#' from 0 to 1. Read each region against its own scale. As with
#' [plot_ordi_star()], the area of the central polygon depends on the order of
#' the spokes and is not an interpretable quantity.
#'
#' @param physeq A phyloseq object containing your microbiome data.
#' @param sample_var A character string. The name of the column in your sample data
#'   to use for grouping samples (e.g., "location", "treatment").
#' @param taxa_rank A character string. The taxonomic rank for the iris plot (default "Phylum").
#' @param n_taxa Integer. The number of top taxa to explicitly show in the iris plot (default 5).
#'   Remaining taxa are grouped into "Other".
#' @param colors_star Optional vector of colors for the central ordination star polygons.
#' @param colors_taxa Optional vector of colors for the outer taxonomic rings.
#' @param view_type Character. Either "separate" (default) or "together".
#'   "separate" facets by group. "together" draws every group into the same
#'   radial band, which overplots and is only readable for a single group.
#' @param error_bar A character string, one of "IQR" (default), "SE", or "none"
#'   for the central star. "IQR" draws the interquartile range of the samples in
#'   the group, around the group **median**; "SE" draws one standard error of
#'   the group **mean**, around the mean.
#' @param fill_alpha Numeric. Transparency of the central star polygon fill.
#' @param distance Ordination distance metric (required if ord is NULL, e.g., "bray").
#' @param ord Optional existing ordination object (e.g., from \code{phyloseq::ordinate()}).
#' @param method A character string specifying the ordination method (default "PCoA").
#'   Constrained methods must be fitted separately and passed via `ord`.
#' @param n_axes Number of axes to use for the central star (default 5).
#' @param plot_order Optional character vector for custom ordering of the sample variable.
#' @param verbose Logical. If TRUE (default), report progress and the
#'   conventions used for the axis percentages.
#'
#' @return A ggplot object representing the star-and-iris visualization.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @import vegan
#' @importFrom stats as.formula sd quantile median
#' @export
#'
#' @examples
#' # data(GlobalPatterns)
#' # gp <- phyloseq::subset_samples(
#' #   GlobalPatterns, SampleType %in% c("Feces", "Ocean", "Soil")
#' # )
#' # plot_ordi_star_iris(
#' #   physeq = gp,
#' #   sample_var = "SampleType",
#' #   taxa_rank = "Phylum",
#' #   distance = "bray"
#' # )
plot_ordi_star_iris <- function(physeq, sample_var, taxa_rank = "Phylum", n_taxa = 5,
                                colors_star = NULL, colors_taxa = NULL,
                                method = "PCoA", view_type = "separate", error_bar = "IQR",
                                fill_alpha = 0.2, distance = NULL, ord = NULL,
                                plot_order = NULL, n_axes = 5, verbose = TRUE) {

  # --- 1. Input Validation ---
  if (!inherits(physeq, "phyloseq")) {
    stop("'physeq' must be a VALID phyloseq object.", call. = FALSE)
  }
  if (!sample_var %in% phyloseq::sample_variables(physeq)) {
    stop(
      "'", sample_var, "' is not a sample variable in the phyloseq object. Available: ",
      paste(phyloseq::sample_variables(physeq), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (taxa_rank != "OTU" && !taxa_rank %in% phyloseq::rank_names(physeq)) {
    stop(
      "'", taxa_rank, "' is not a taxonomic rank in the phyloseq object. Available: ",
      paste(phyloseq::rank_names(physeq), collapse = ", "), ".",
      call. = FALSE
    )
  }
  # Check for ggnewscale
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("The 'ggnewscale' package is required for this beta function. Please install it using install.packages('ggnewscale').", call. = FALSE)
  }

  view_type <- match.arg(view_type, c("together", "separate"))
  error_bar <- match.arg(error_bar, c("IQR", "SE", "none"))

  # --- 2. Ordination ---
  method <- match.arg(method, ordi_methods())

  if (is.null(ord)) {
    if (method %in% c("CCA", "RDA", "CAP")) {
      stop(
        method, " is a constrained ordination and needs a formula, which cannot ",
        "be supplied through this function. Fit it with phyloseq::ordinate() ",
        "and pass the result via the 'ord' argument.",
        call. = FALSE
      )
    }
    if (is.null(distance)) {
      stop("'distance' must be provided for ", method, " when 'ord' is NULL.", call. = FALSE)
    }
    if (verbose) message("Calculating dissimilarity and performing ", method, " for the central star...")
    ord <- phyloseq::ordinate(physeq, method = method, distance = distance)
  } else if (verbose) {
    message("Using provided ordination object...")
  }

  axes <- extract_ordination_axes(ord, n_axes = n_axes, verbose = verbose)

  if (axes$n_axes < n_axes) {
    warning(method, " produced fewer than ", n_axes, " axes. Plotting the ",
            axes$n_axes, " available.", call. = FALSE)
  }
  n_axes <- axes$n_axes
  axis_labels <- axes$labels

  # Labels are assigned after conversion: they contain characters make.names()
  # would mangle.
  pcoa_axes <- as.data.frame(axes$scores)
  colnames(pcoa_axes) <- axis_labels

  pcoa_df <- cbind(pcoa_axes, align_sample_data(axes$scores, physeq))

  y_offset <- radial_offset(pcoa_axes)
  pcoa_df_shifted <- pcoa_df
  pcoa_df_shifted[axis_labels] <- pcoa_df[axis_labels] + y_offset

  pcoa_stats <- pcoa_df_shifted %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(axis_labels), names_to = "Axis", values_to = "Value") %>%
    dplyr::group_by(!!dplyr::sym(sample_var), Axis) %>%
    dplyr::summarise(
      N_samples = sum(!is.na(Value)),
      Mean_Position = mean(Value, na.rm = TRUE),
      Median_Position = stats::median(Value, na.rm = TRUE),
      Q25_Position = stats::quantile(Value, 0.25, na.rm = TRUE),
      Q75_Position = stats::quantile(Value, 0.75, na.rm = TRUE),
      SE_Position = stats::sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SE_Min = Mean_Position - SE_Position,
      SE_Max = Mean_Position + SE_Position
    )

  # The IQR brackets the median; the SE describes the mean. Plot whichever
  # centre matches the interval being drawn.
  center_col <- if (error_bar == "IQR") "Median_Position" else "Mean_Position"
  pcoa_stats$Center_Position <- pcoa_stats[[center_col]]
  center_label <- if (error_bar == "IQR") "Median" else "Mean"

  if (error_bar != "none") {
    thin <- unique(as.character(pcoa_stats[[sample_var]][pcoa_stats$N_samples < 2]))
    if (length(thin) > 0) {
      warning(
        "No interval can be drawn for group(s) with a single sample: ",
        paste(thin, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  pcoa_stats$Axis <- factor(pcoa_stats$Axis, levels = axis_labels)
  pcoa_stats$Axis_numeric <- as.numeric(pcoa_stats$Axis)

  # --- 3. Taxonomic Aggregation (Iris Plot Data) ---
  if (verbose) message("Aggregating taxonomy to ", taxa_rank, " level...")
  physeq_rel <- phyloseq::transform_sample_counts(physeq, to_relative)

  if (taxa_rank != "OTU") {
    physeq_rel <- phyloseq::tax_glom(physeq_rel, taxrank = taxa_rank)
  }

  df_tax <- phyloseq::psmelt(physeq_rel)

  # Determine Top Taxa overall in the object
  top_taxa <- df_tax %>%
    dplyr::group_by(!!dplyr::sym(taxa_rank)) %>%
    dplyr::summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(TotalAbundance)) %>%
    dplyr::slice_head(n = n_taxa) %>%
    dplyr::pull(!!dplyr::sym(taxa_rank)) %>%
    as.character()

  df_tax <- df_tax %>%
    dplyr::mutate(Taxa_Group = dplyr::if_else(as.character(get(taxa_rank)) %in% top_taxa, as.character(get(taxa_rank)), "Other"))

  # Summarize by sample and new Taxa_Group
  df_tax_samples <- df_tax %>%
    dplyr::group_by(Sample, !!dplyr::sym(sample_var), Taxa_Group) %>%
    dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

  # Stack the named taxa in abundance order with the pooled wedge outermost,
  # rather than in whatever order the labels happen to sort into.
  df_tax_samples$Taxa_Group <- order_taxa_group(df_tax_samples$Taxa_Group, levels_order = top_taxa)

  # Calculate mapping parameters for the ring.
  total_circumference <- n_axes

  df_tax_pos <- df_tax_samples %>%
    dplyr::distinct(!!dplyr::sym(sample_var), Sample) %>%
    dplyr::group_by(!!dplyr::sym(sample_var)) %>%
    dplyr::arrange(Sample, .by_group = TRUE) %>%
    dplyr::mutate(
      N_ring_samples = dplyr::n(),
      sample_idx = dplyr::row_number(),
      width = total_circumference / N_ring_samples,
      x_pos = 0.5 + width * (sample_idx - 0.5)
    ) %>%
    dplyr::ungroup()

  # Calculate radii (Y-axis scaling). Use only the elements actually drawn, so
  # error_bar = "none" does not leave a gap sized for an absent interval.
  star_cols <- switch(error_bar,
    "IQR" = c("Center_Position", "Q75_Position"),
    "SE" = c("Center_Position", "SE_Max"),
    "none" = "Center_Position"
  )
  max_star_y <- suppressWarnings(max(unlist(pcoa_stats[star_cols]), na.rm = TRUE))
  if (!is.finite(max_star_y)) max_star_y <- max(pcoa_stats$Center_Position, na.rm = TRUE)

  # Empty center ends slightly above the central star elements
  r_inner <- max_star_y * 1.3
  r_thickness <- max_star_y * 0.8

  iris_df <- df_tax_samples %>%
    dplyr::left_join(df_tax_pos %>% dplyr::select(Sample, x_pos, width), by = "Sample") %>%
    dplyr::group_by(!!dplyr::sym(sample_var), Sample) %>%
    dplyr::arrange(Taxa_Group, .by_group = TRUE) %>%
    dplyr::mutate(
      cum_abundance = cumsum(Abundance),
      ymax = r_inner + cum_abundance * r_thickness,
      ymin = r_inner + (cum_abundance - Abundance) * r_thickness,
      xmin = x_pos - width / 2,
      xmax = x_pos + width / 2
    ) %>%
    dplyr::ungroup()

  # Convert geom_rect data into explicit polygons to avoid coord_radar()
  # conversion issues. Built vectorised: one row per corner, four corners per
  # rectangle, rather than one data.frame per rectangle.
  n_rect <- nrow(iris_df)
  corner <- rep(seq_len(n_rect), each = 4L)
  iris_poly_df <- data.frame(
    Sample = iris_df$Sample[corner],
    Taxa_Group = iris_df$Taxa_Group[corner],
    poly_x = as.vector(t(cbind(iris_df$xmin, iris_df$xmin, iris_df$xmax, iris_df$xmax))),
    poly_y = as.vector(t(cbind(iris_df$ymin, iris_df$ymax, iris_df$ymax, iris_df$ymin))),
    poly_group = paste0(iris_df$Sample[corner], "_", iris_df$Taxa_Group[corner], "_", corner),
    stringsAsFactors = FALSE
  )
  iris_poly_df[[sample_var]] <- iris_df[[sample_var]][corner]

  # --- 4. Plot Ordering ---
  if (!is.null(plot_order)) {
    if (length(plot_order) == 1 && identical(plot_order, "hclust")) {
      stop(
        "plot_order = \"hclust\" is not supported by plot_ordi_star_iris(); ",
        "supply an explicit ordering vector instead.",
        call. = FALSE
      )
    }
    pcoa_stats[[sample_var]] <- factor(pcoa_stats[[sample_var]], levels = plot_order)
    iris_poly_df[[sample_var]] <- factor(iris_poly_df[[sample_var]], levels = plot_order)
  }

  # --- 5. Colors ---
  groups <- unique(as.character(pcoa_stats[[sample_var]]))
  if (is.null(colors_star)) {
    colors_star <- get_default_colors(groups)
  }

  taxa_groups_all <- levels(droplevels(as.factor(iris_poly_df$Taxa_Group)))
  if (is.null(colors_taxa)) {
    colors_taxa <- get_default_colors(taxa_groups_all)
  }

  # --- 6. Plot Construction ---
  title_suffix <- switch(error_bar, "IQR" = " with IQR", "SE" = " with SE", "none" = "")
  caption_text <- switch(error_bar,
    "IQR" = "Centre: bars span the interquartile range of samples within each group, around the group median. Outer ring: relative abundance, 0 to 1, one segment per sample.",
    "SE" = "Centre: bars span one standard error of the group mean. Outer ring: relative abundance, 0 to 1, one segment per sample.",
    "none" = "Outer ring: relative abundance, 0 to 1, one segment per sample."
  )

  p <- ggplot2::ggplot() +
    # Layer 1: Ordination Star Polygon
    ggplot2::geom_polygon(data = pcoa_stats, ggplot2::aes(x = Axis_numeric, y = Center_Position, group = !!dplyr::sym(sample_var), fill = !!dplyr::sym(sample_var), color = !!dplyr::sym(sample_var)), linewidth = 0.8, alpha = fill_alpha)

  # Error bars for the central star
  if (error_bar == "IQR") {
    p <- p + ggplot2::geom_errorbar(data = pcoa_stats, ggplot2::aes(x = Axis_numeric, ymin = Q25_Position, ymax = Q75_Position, color = !!dplyr::sym(sample_var)), width = 0.2, alpha = 0.7, show.legend = FALSE)
  } else if (error_bar == "SE") {
    p <- p + ggplot2::geom_errorbar(data = pcoa_stats, ggplot2::aes(x = Axis_numeric, ymin = SE_Min, ymax = SE_Max, color = !!dplyr::sym(sample_var)), width = 0.2, alpha = 0.7, show.legend = FALSE)
  }

  # Set scales for Layer 1
  p <- p +
    ggplot2::scale_fill_manual(values = colors_star, name = paste(method, "Star\n", sample_var)) +
    ggplot2::scale_color_manual(values = colors_star, name = paste(method, "Star\n", sample_var)) +

    # NEW FILL SCALE for the Outer Ring
    ggnewscale::new_scale_fill() +

    # Layer 2: Taxonomic Ring (Iris) mapped via geom_polygon
    ggplot2::geom_polygon(data = iris_poly_df, ggplot2::aes(x = poly_x, y = poly_y, group = poly_group, fill = Taxa_Group), color = "white", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = colors_taxa, name = paste("Outer Iris\n", taxa_rank)) +

    # Theme & Coordinates
    ggplot2::scale_x_continuous(breaks = 1:n_axes, labels = axis_labels, limits = c(0.5, n_axes + 0.5)) +
    # Radial labels apply to the central star only; the ring has its own 0-1
    # scale, so labelling it with ordination scores would be misleading.
    ggplot2::scale_y_continuous(
      labels = function(x) ifelse(x > r_inner, "", sprintf("%.2f", x - y_offset)),
      breaks = scales::breaks_pretty(n = 3)
    ) +
    coord_radar(clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "serif"),
      axis.text.x = ggplot2::element_text(color = "black", size = 10),
      axis.text.y = ggplot2::element_text(color = "gray50"),
      panel.grid.major = ggplot2::element_line(color = "#e8e8e8", linewidth = 0.5),
      plot.caption = ggplot2::element_text(color = "gray40", hjust = 0, size = 8),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    ) +
    ggplot2::labs(
      title = paste0(center_label, " ", method, " Position", title_suffix, " & Taxa Iris"),
      x = "",
      y = paste0(center_label, " Score (centre)"),
      caption = caption_text
    )

  if (view_type == "separate") {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", sample_var)))
  } else {
    warning(
      "view_type = \"together\" draws every group into the same radial band, ",
      "which overplots when there is more than one group.",
      call. = FALSE
    )
  }

  return(p)
}
