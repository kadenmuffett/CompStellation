#' Plot a Star Plot of Mean PCoA Position for Groups
#'
#' This function calculates dissimilarity, performs PCoA, and plots the
#' mean position of sample groups along the first 5 PCoA axes using a star (radar) plot.
#'
#' @param physeq A phyloseq object containing your microbiome data.
#' @param sample_var A character string. The name of the column in your sample data
#'   to use for grouping samples (e.g., "location", "treatment").
#' @param view_type A character string, either "together" (default) or "separate".
#'   Determines whether groups are plotted on the same plot or faceted.
#' @param error_bar A character string, one of "IQR" (default), "SE", or "none".
#'   Determines the type of error bars to display.
#' @param fill_alpha A numeric value between 0 and 1 (default 0.2).
#'   Controls the transparency of the polygon fill under the star plot.
#' @param distance An accepted phyloseq ordination type ("bray" for example). Required if `ord` is NULL.
#' @param ord An optional existing ordination object (e.g., from \code{phyloseq::ordinate()}). If provided, ordination is skipped.
#' @param plot_order A character vector for custom ordering of the sample variable, "hclust" for Ward's clustering based on Euclidean distance, or NULL (default) for alphabetical.
#'
#' @return A ggplot object representing the star plot of PCoA centroids.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @import vegan
#'
#' @importFrom stats as.formula sd quantile
#'
#' @export
#'
#' @examples
#' # ...
#' #   plot_ordi_star(
#' #     physeq = gp_subset,
#' #     sample_var = "SampleType",
#' #     colors_all = c(...),
#' #     view_type = "together",
#' #     error_bar = "SE",
#' #     fill_alpha = 0.2,
#' #     distance = "bray"
#' #    )
#' #
plot_ordi_star <- function(physeq, sample_var, colors_all, method = "PCoA", view_type = "together", error_bar = "IQR", fill_alpha = 0.2, distance = NULL, ord = NULL, plot_order = NULL, n_axes = 5) {
  # --- 1. Input Validation ---
  if (!inherits(physeq, "phyloseq")) {
    stop("Error: 'physeq' must be a VALID phyloseq object.")
  }
  if (!sample_var %in% phyloseq::sample_variables(physeq)) {
    stop("Error: '", sample_var, "' is not a valid sample variable in the phyloseq object.")
  }

  view_type <- match.arg(view_type, c("together", "separate"))
  error_bar <- match.arg(error_bar, c("IQR", "SE", "none"))

  # --- 2. Ordination ---
  method <- match.arg(method, c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"))

  if (is.null(ord)) {
    if (is.null(distance) && method != "PCA") {
      stop("Error: 'distance' must be provided for ", method, " if 'ord' is NULL.")
    }
    
    message("Calculating dissimilarity and performing ", method, "...")
    if (method == "PCA") {
      ord <- phyloseq::ordinate(physeq, method = "PCA")
    } else {
      ord <- phyloseq::ordinate(physeq, method = method, distance = distance)
    }
  } else {
    message("Using provided ordination object...")
  }

  # Generic extraction logic
  if (inherits(ord, c("pcoa", "dpcoa", "list")) && !is.null(ord$vectors)) {
    vecs <- ord$vectors
    eigen_vals <- ord$values$Eigenvalues
    if(!is.null(eigen_vals)) percent_var <- 100 * (eigen_vals / sum(eigen_vals)) else percent_var <- NULL
    axis_raw_names <- paste0("Axis ", 1:ncol(vecs))
  } else if (inherits(ord, c("metaMDS", "monoMDS"))) {
    vecs <- ord$points
    percent_var <- NULL
    axis_raw_names <- paste0("NMDS", 1:ncol(vecs))
  } else if (inherits(ord, c("rda", "cca", "prcomp"))) {
    vecs <- vegan::scores(ord, display = "sites")
    eigen_vals <- ord$CA$eig
    if(!is.null(eigen_vals)) percent_var <- 100 * (eigen_vals / sum(eigen_vals)) else percent_var <- NULL
    axis_raw_names <- paste0("PC", 1:ncol(vecs))
  } else {
    vecs <- tryCatch(vegan::scores(ord, display = "sites"), error = function(e) ord$vectors)
    percent_var <- NULL
    axis_raw_names <- paste0("Axis ", 1:ncol(vecs))
  }

  actual_axes <- min(n_axes, ncol(vecs))
  if (actual_axes < n_axes) {
    warning(method, " resulted in fewer than ", n_axes, " axes. Plotting ", actual_axes, " available axes.")
  }
  n_axes <- actual_axes

  if (!is.null(percent_var)) {
    axis_labels <- paste0(axis_raw_names[1:n_axes], "\n(", sprintf("%.1f", percent_var[1:n_axes]), "%)")
  } else {
    axis_labels <- axis_raw_names[1:n_axes]
  }

  pcoa_axes <- data.frame(vecs[, 1:n_axes, drop = FALSE])
  colnames(pcoa_axes) <- axis_labels # Use the new labels as column names immediately

  # Add sample data
  sample_dat <- data.frame(phyloseq::sample_data(physeq))
  pcoa_df <- cbind(pcoa_axes, sample_dat)

  # --- 3. Aggregation (Calculate Mean Position and stats) ---

  # Handle Negative Values for Radar Plot
  # Shift all values to be positive if necessary
  # Find global minimum across all axes
  min_val <- min(pcoa_axes, na.rm = TRUE)

  if (min_val < 0) {
    # Add a buffer (e.g., 5% of range or just absolute min)
    # Let's just shift so min becomes 0 (or slightly above 0 to avoid touching center?)
    # Center of radar is 0. If data is 0, it vanishes.
    # Let's shift so min is e.g. 0.1? Or just 0.
    y_offset <- abs(min_val) * 1.1 # Shift by 110% of min to ensure it's clearly positive
  } else {
    y_offset <- 0
  }

  # Shift the axis columns
  pcoa_df_shifted <- pcoa_df
  pcoa_df_shifted[axis_labels] <- pcoa_df[axis_labels] + y_offset

  # Group by sample_var and calculate stats for each axis
  # We use the new column names (axis_labels)

  # This is tricky with 'across' returning multiple values.
  # Better approach: Pivot longer first, then summarise.
  pcoa_long_raw <- pcoa_df_shifted %>%
    tidyr::pivot_longer(cols = all_of(axis_labels), names_to = "Axis", values_to = "Value")

  pcoa_stats <- pcoa_long_raw %>%
    group_by(!!sym(sample_var), Axis) %>%
    summarise(
      Mean_Position = mean(Value, na.rm = TRUE),
      # IQR Stats
      Q25_Position = quantile(Value, 0.25, na.rm = TRUE),
      Q75_Position = quantile(Value, 0.75, na.rm = TRUE),
      # SE Stats
      SE_Position = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      SE_Min = Mean_Position - SE_Position,
      SE_Max = Mean_Position + SE_Position
    )

  # Ensure Axis order is preserved
  pcoa_stats$Axis <- factor(pcoa_stats$Axis, levels = axis_labels)

  # --- 3b. Plot Ordering ---
  if (!is.null(plot_order)) {
    if (length(plot_order) == 1 && plot_order == "hclust") {
      # Calculate hclust
      wide_df <- pcoa_stats %>%
        dplyr::select(!!sym(sample_var), Axis, Mean_Position) %>%
        tidyr::pivot_wider(names_from = Axis, values_from = Mean_Position, values_fill = list(Mean_Position = 0))

      dist_mat <- vegan::vegdist(wide_df %>% dplyr::select(-!!sym(sample_var)), method = "bray")
      hc <- stats::hclust(dist_mat, method = "complete")
      ordered_groups <- wide_df[[sample_var]][hc$order]

      pcoa_stats[[sample_var]] <- factor(pcoa_stats[[sample_var]], levels = ordered_groups)
    } else {
      # Custom order provided by user
      pcoa_stats[[sample_var]] <- factor(pcoa_stats[[sample_var]], levels = plot_order)
    }
  }

  # --- 4. Plotting ---

  # Handle default colors if missing
  if (missing(colors_all) || (length(colors_all) == 1 && colors_all == "hclust")) {
    if (is.factor(pcoa_stats[[sample_var]])) {
      groups <- levels(pcoa_stats[[sample_var]])
    } else {
      groups <- as.character(unique(pcoa_stats[[sample_var]]))
    }

    if (!missing(colors_all) && colors_all == "hclust") {
      long_df <- pcoa_stats %>%
        dplyr::select(!!sym(sample_var), Axis, Mean_Position) %>%
        tidyr::pivot_wider(names_from = Axis, values_from = Mean_Position, values_fill = list(Mean_Position = 0))

      dist_m <- vegan::vegdist(long_df %>% dplyr::select(-!!sym(sample_var)), method = "bray")
      hc_res <- stats::hclust(dist_m, method = "complete")
      ordered_names <- long_df[[sample_var]][hc_res$order]

      colors_all <- get_hclust_colors(hc_res, ordered_names)
    } else {
      colors_all <- get_default_colors(groups)
    }
  }

  title_suffix <- switch(error_bar,
    "IQR" = " with IQR",
    "SE" = " with SE",
    "none" = ""
  )

  star_plot <- ggplot(pcoa_stats, aes(
    x = Axis, y = Mean_Position, group = !!sym(sample_var),
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
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(
      title = paste0("Mean ", method, " Position (First ", n_axes, " Axes)", title_suffix),
      x = "",
      y = "Mean Score",
      color = sample_var,
      fill = sample_var
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
