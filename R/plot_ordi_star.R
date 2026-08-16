#' Plot a Star Plot of Group Position Across Several Ordination Axes
#'
#' Runs an ordination (PCoA by default) on a beta-diversity distance matrix and
#' summarises each group's position by plotting its centre on the first
#' `n_axes` axes as the spokes of a star (radar) plot. This lets several axes be
#' compared at once, rather than the two a conventional ordination scatterplot
#' shows.
#'
#' Individual samples are not drawn, and the plot carries no taxonomic
#' information. For ordination position alongside taxonomic composition see
#' [plot_ordi_star_iris()].
#'
#' @section Interpreting the plot:
#' The radial axis is **offset**. Ordination scores are negative as well as
#' positive, so every axis is shifted by a single constant to make the values
#' plottable on a radius. Radii are therefore not proportional to the scores and
#' are not comparable between axes; the tick labels show the original,
#' unshifted values. The area enclosed by a polygon depends on the order of the
#' spokes and is not an interpretable quantity. The plot is descriptive: it does
#' not test whether groups differ, so pair it with an appropriate test such as
#' [vegan::adonis2()] on the same distance matrix.
#'
#' @param physeq A phyloseq object containing your microbiome data.
#' @param sample_var A character string. The name of the column in your sample data
#'   to use for grouping samples (e.g., "location", "treatment").
#' @param colors_all Optional. A character vector of colors, or "hclust" for automatic clustering colors.
#' @param method A character string specifying the ordination method (default "PCoA").
#'   One of "PCoA", "MDS", "NMDS", "DPCoA", "DCA", "CCA", "RDA" or "CAP".
#'   Constrained methods ("CCA", "RDA", "CAP") need a formula, so they must be
#'   fitted separately and passed via `ord`.
#' @param view_type A character string, either "together" (default) or "separate".
#'   Determines whether groups are plotted on the same plot or faceted.
#' @param error_bar A character string, one of "IQR" (default), "SE", or "none".
#'   "IQR" draws the interquartile range of the samples in the group, around the
#'   group **median**; "SE" draws one standard error of the group **mean**,
#'   around the mean. The two describe different things — spread of the samples
#'   versus precision of the centre — and the plot title records which was used.
#' @param fill_alpha A numeric value between 0 and 1 (default 0.2).
#'   Controls the transparency of the polygon fill under the star plot.
#' @param distance An accepted phyloseq ordination type ("bray" for example). Required if `ord` is NULL.
#' @param ord An optional existing ordination object (e.g., from \code{phyloseq::ordinate()}). If provided, ordination is skipped.
#' @param plot_order A character vector for custom ordering of the sample variable,
#'   "hclust" for complete-linkage clustering on Euclidean distances between the
#'   group centres, or NULL (default) for alphabetical.
#' @param n_axes Number of axes to plot.
#' @param base_colors Optional. A character vector of base colors to use for hclust coloring.
#' @param verbose Logical. If TRUE (default), report progress and the
#'   conventions used for the axis percentages.
#'
#' @return A ggplot object representing the star plot of group ordination centres.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @import vegan
#'
#' @importFrom stats as.formula sd quantile median dist
#'
#' @export
#'
#' @examples
#' # data(GlobalPatterns)
#' # gp <- phyloseq::subset_samples(
#' #   GlobalPatterns, SampleType %in% c("Feces", "Ocean", "Soil")
#' # )
#' # plot_ordi_star(
#' #   physeq = gp,
#' #   sample_var = "SampleType",
#' #   distance = "bray",
#' #   view_type = "together",
#' #   error_bar = "SE"
#' # )
plot_ordi_star <- function(physeq, sample_var, colors_all, method = "PCoA", view_type = "together", error_bar = "IQR", fill_alpha = 0.2, distance = NULL, ord = NULL, plot_order = NULL, n_axes = 5, base_colors = NULL, verbose = TRUE) {
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

    if (verbose) message("Calculating dissimilarity and performing ", method, "...")
    ord <- phyloseq::ordinate(physeq, method = method, distance = distance)
  } else if (verbose) {
    message("Using provided ordination object...")
  }

  axes <- extract_ordination_axes(ord, n_axes = n_axes, verbose = verbose)

  if (axes$n_axes < n_axes) {
    warning(
      method, " produced fewer than ", n_axes, " axes. Plotting the ",
      axes$n_axes, " available.",
      call. = FALSE
    )
  }
  n_axes <- axes$n_axes
  axis_labels <- axes$labels

  # Assign the labels after the data.frame conversion: they contain newlines and
  # percent signs that make.names() would mangle.
  pcoa_axes <- as.data.frame(axes$scores)
  colnames(pcoa_axes) <- axis_labels

  # Match sample data on sample names rather than binding by position
  pcoa_df <- cbind(pcoa_axes, align_sample_data(axes$scores, physeq))

  # --- 3. Aggregation (Calculate Group Centres and Spread) ---

  # A radar radius cannot show negative values, so shift every axis by one
  # constant. Because the shift is common to all axes it cancels out of any
  # Euclidean distance computed below, and the axis labels undo it for display.
  y_offset <- radial_offset(pcoa_axes)

  pcoa_df_shifted <- pcoa_df
  pcoa_df_shifted[axis_labels] <- pcoa_df[axis_labels] + y_offset

  pcoa_long_raw <- pcoa_df_shifted %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(axis_labels), names_to = "Axis", values_to = "Value")

  pcoa_stats <- pcoa_long_raw %>%
    group_by(!!sym(sample_var), Axis) %>%
    summarise(
      N_samples = sum(!is.na(Value)),
      Mean_Position = mean(Value, na.rm = TRUE),
      Median_Position = stats::median(Value, na.rm = TRUE),
      # IQR Stats
      Q25_Position = quantile(Value, 0.25, na.rm = TRUE),
      Q75_Position = quantile(Value, 0.75, na.rm = TRUE),
      # SE Stats
      SE_Position = sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
      .groups = "drop"
    ) %>%
    mutate(
      SE_Min = Mean_Position - SE_Position,
      SE_Max = Mean_Position + SE_Position
    )

  # The IQR brackets the median, so plot the median when the IQR is shown;
  # the SE describes the mean, so plot the mean when the SE is shown.
  center_col <- if (error_bar == "IQR") "Median_Position" else "Mean_Position"
  pcoa_stats$Center_Position <- pcoa_stats[[center_col]]
  center_label <- if (error_bar == "IQR") "Median" else "Mean"

  # Groups of one have no spread to show; say so rather than dropping the bar
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

  # Ensure Axis order is preserved. The rows must be sorted as well as the
  # levels set: geom_polygon() joins vertices in row order, so a factor whose
  # levels disagree with the row order produces a self-crossing shape. The rows
  # arrive sorted alphabetically by label, which stops matching the axis numbers
  # at ten or more axes ("Axis 10" sorts before "Axis 2").
  pcoa_stats$Axis <- factor(pcoa_stats$Axis, levels = axis_labels)
  pcoa_stats <- pcoa_stats %>% dplyr::arrange(!!sym(sample_var), Axis)

  # --- 3b. Plot Ordering ---
  # Ordination scores are Euclidean by construction, and are signed. Bray-Curtis
  # presupposes non-negative abundance-like data, so it is not a valid distance
  # here; Euclidean distance is, and it is unaffected by the radial offset.
  group_centres <- function() {
    pcoa_stats %>%
      dplyr::select(!!sym(sample_var), Axis, Center_Position) %>%
      tidyr::pivot_wider(names_from = Axis, values_from = Center_Position, values_fill = list(Center_Position = 0))
  }

  cluster_groups <- function(wide_df) {
    dist_mat <- stats::dist(
      as.matrix(wide_df %>% dplyr::select(-!!sym(sample_var))),
      method = "euclidean"
    )
    hc <- stats::hclust(dist_mat, method = "complete")
    hc$labels <- as.character(wide_df[[sample_var]])
    hc
  }

  if (!is.null(plot_order)) {
    if (length(plot_order) == 1 && plot_order == "hclust") {
      if (verbose) {
        message("Ordering groups by complete-linkage clustering on Euclidean distances between group centres...")
      }
      wide_df <- group_centres()
      hc <- cluster_groups(wide_df)
      ordered_groups <- hc$labels[hc$order]

      pcoa_stats[[sample_var]] <- factor(pcoa_stats[[sample_var]], levels = ordered_groups)
    } else {
      # Custom order provided by user
      pcoa_stats[[sample_var]] <- factor(pcoa_stats[[sample_var]], levels = plot_order)
    }
  }

  # --- 4. Plotting ---

  # Handle default colors if missing
  if (missing(colors_all) || is.null(colors_all) || (length(colors_all) == 1 && colors_all == "hclust")) {
    if (is.factor(pcoa_stats[[sample_var]])) {
      groups <- levels(pcoa_stats[[sample_var]])
    } else {
      groups <- as.character(unique(pcoa_stats[[sample_var]]))
    }

    if (!missing(colors_all) && !is.null(colors_all) && colors_all == "hclust") {
      hc_res <- cluster_groups(group_centres())
      ordered_names <- hc_res$labels[hc_res$order]

      colors_all <- get_hclust_colors(hc_res, ordered_names, base_colors = base_colors)
    } else {
      colors_all <- get_default_colors(groups)
    }
  } else if (!is.null(plot_order) && length(plot_order) == 1 && plot_order == "hclust") {
    # If hclust + color options specified, treat colors as bases for hclust base_palette
    hc_res <- cluster_groups(group_centres())
    ordered_names <- hc_res$labels[hc_res$order]

    colors_all <- get_hclust_colors(hc_res, ordered_names, base_colors = colors_all)
  }

  title_suffix <- switch(error_bar,
    "IQR" = " with IQR",
    "SE" = " with SE",
    "none" = ""
  )

  caption_text <- switch(error_bar,
    "IQR" = "Bars span the interquartile range of samples within each group, around the group median.",
    "SE" = "Bars span one standard error of the group mean.",
    "none" = NULL
  )

  star_plot <- ggplot(pcoa_stats, aes(
    x = Axis, y = Center_Position, group = !!sym(sample_var),
    color = !!sym(sample_var), fill = !!sym(sample_var)
  )) +
    geom_polygon(aes(), linewidth = 0.8, alpha = fill_alpha) +
    scale_fill_manual(values = colors_all) +
    scale_color_manual(values = colors_all) +
    # Adjust Y-axis labels to show original values
    scale_y_continuous(labels = function(x) sprintf("%.2f", x - y_offset), breaks = scales::breaks_pretty(n = 5)) +
    theme_minimal() +
    theme(
      text = element_text(family = "serif"), axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "gray50"),
      panel.grid.major = element_line(color = "#e8e8e8", linewidth = 0.5),
      plot.caption = element_text(color = "gray40", hjust = 0, size = 8),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(
      title = paste0(center_label, " ", method, " Position (First ", n_axes, " Axes)", title_suffix),
      x = "",
      y = paste0(center_label, " Score"),
      color = sample_var,
      fill = sample_var,
      caption = caption_text
    ) +
    coord_radar(clip = "off") # Uses the function from utils-radar.R

  # Add Error Bars conditionally
  if (error_bar == "IQR") {
    star_plot <- star_plot +
      geom_errorbar(aes(ymin = Q25_Position, ymax = Q75_Position), width = 0.2, alpha = 0.7)
  } else if (error_bar == "SE") {
    star_plot <- star_plot +
      geom_errorbar(aes(ymin = SE_Min, ymax = SE_Max), width = 0.2, alpha = 0.7)
  }

  if (view_type == "separate") {
    star_plot <- star_plot + facet_wrap(as.formula(paste("~", sample_var)))
  }

  return(star_plot)
}
