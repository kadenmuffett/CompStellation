#' Plot a Star Plot of Mean PCoA Position Surrounded by a Taxonomic Iris Plot
#'
#' This experimental (beta) function combines a 5-axis PCoA radar plot (the "star")
#' in the center, and surrounds it with a circular stacked bar chart (the "iris plot")
#' representing the taxonomic composition of all individual samples belonging to that group.
#'
#' @param physeq A phyloseq object containing your microbiome data.
#' @param sample_var A character string. The name of the column in your sample data
#'   to use for grouping samples (e.g., "location", "treatment").
#' @param taxa_rank A character string. The taxonomic rank for the iris plot (default "Phylum").
#' @param n_taxa Integer. The number of top taxa to explicitly show in the iris plot (default 5).
#'   Remaining taxa are grouped into "Other".
#' @param colors_star Optional vector of colors for the central PCoA star polygons.
#' @param colors_taxa Optional vector of colors for the outer taxonomic rings.
#' @param view_type Character. Either "separate" (default) or "together".
#'   Note: "together" concentrically stacks rings and polygons, which may be unreadable
#'   with many groups or many taxa. "separate" facets by group for clarity.
#' @param error_bar A character string, one of "IQR" (default), "SE", or "none" for the PCoA star.
#' @param fill_alpha Numeric. Transparency of the central star polygon fill.
#' @param distance Ordination distance metric (required if ord is NULL, e.g., "bray").
#' @param ord Optional existing ordination object (e.g., from \code{phyloseq::ordinate()}).
#' @param plot_order Optional character vector for custom ordering of the sample variable.
#'
#' @return A ggplot object representing the scatter-iris visualization.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @import vegan
#' @importFrom stats as.formula sd quantile
#' @export
plot_ordi_star_iris <- function(physeq, sample_var, taxa_rank = "Phylum", n_taxa = 5,
                                colors_star = NULL, colors_taxa = NULL,
                                method = "PCoA", view_type = "separate", error_bar = "IQR",
                                fill_alpha = 0.2, distance = NULL, ord = NULL,
                                plot_order = NULL, n_axes = 5) {

  # --- 1. Input Validation ---
  if (!inherits(physeq, "phyloseq")) {
    stop("Error: 'physeq' must be a VALID phyloseq object.")
  }
  if (!sample_var %in% phyloseq::sample_variables(physeq)) {
    stop("Error: '", sample_var, "' is not a valid sample variable in the phyloseq object.")
  }
  # Check for ggnewscale
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("The 'ggnewscale' package is required for this beta function. Please install it using install.packages('ggnewscale').")
  }

  view_type <- match.arg(view_type, c("together", "separate"))
  error_bar <- match.arg(error_bar, c("IQR", "SE", "none"))

  # --- 2. PCoA Calculation ---
  method <- match.arg(method, c("PCoA", "NMDS", "PCA"))
  
  if (is.null(ord)) {
    if (is.null(distance) && method != "PCA") stop("Error: 'distance' must be provided for ", method, " if 'ord' is NULL.")
    message("Calculating dissimilarity and performing ", method, " for the central star...")
    if (method == "PCA") {
      ord <- phyloseq::ordinate(physeq, method = "PCA")
    } else {
      ord <- phyloseq::ordinate(physeq, method = method, distance = distance)
    }
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
    warning(method, " resulted in fewer than ", n_axes, " axes.")
  }
  n_axes <- actual_axes

  if (!is.null(percent_var)) {
    axis_labels <- paste0(axis_raw_names[1:n_axes], "\n(", sprintf("%.1f", percent_var[1:n_axes]), "%)")
  } else {
    axis_labels <- axis_raw_names[1:n_axes]
  }

  pcoa_axes <- data.frame(vecs[, 1:n_axes, drop = FALSE])
  colnames(pcoa_axes) <- axis_labels
  sample_dat <- data.frame(phyloseq::sample_data(physeq))
  pcoa_df <- cbind(pcoa_axes, sample_dat)
  pcoa_df$Sample_ID_unique <- phyloseq::sample_names(physeq)

  min_val <- min(pcoa_axes, na.rm = TRUE)
  y_offset <- ifelse(min_val < 0, abs(min_val) * 1.1, 0)
  pcoa_df_shifted <- pcoa_df
  pcoa_df_shifted[axis_labels] <- pcoa_df[axis_labels] + y_offset

  pcoa_stats <- pcoa_df_shifted %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(axis_labels), names_to = "Axis", values_to = "Value") %>%
    dplyr::group_by(!!dplyr::sym(sample_var), Axis) %>%
    dplyr::summarise(
      Mean_Position = mean(Value, na.rm = TRUE),
      Q25_Position = stats::quantile(Value, 0.25, na.rm = TRUE),
      Q75_Position = stats::quantile(Value, 0.75, na.rm = TRUE),
      SE_Position = stats::sd(Value, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SE_Min = Mean_Position - SE_Position,
      SE_Max = Mean_Position + SE_Position
    )

  pcoa_stats$Axis <- factor(pcoa_stats$Axis, levels = axis_labels)
  pcoa_stats$Axis_numeric <- as.numeric(pcoa_stats$Axis)

  # --- 3. Taxonomic Aggregation (Iris Plot Data) ---
  message(paste0("Aggregating taxonomy to ", taxa_rank, " level..."))
  physeq_rel <- phyloseq::transform_sample_counts(physeq, function(x) 1 * x / sum(x))
  
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

  # Calculate mapping parameters for the ring.
  total_circumference <- n_axes

  df_tax_pos <- df_tax_samples %>%
    dplyr::group_by(!!dplyr::sym(sample_var), Sample) %>%
    dplyr::summarise(TotalAbun = sum(Abundance, na.rm=TRUE), .groups = "drop") %>%
    dplyr::group_by(!!dplyr::sym(sample_var)) %>%
    dplyr::arrange(Sample, .by_group = TRUE) %>%
    dplyr::mutate(
      N_samples = dplyr::n(),
      sample_idx = dplyr::row_number(),
      width = total_circumference / N_samples, 
      x_pos = 0.5 + width * (sample_idx - 0.5)
    ) %>%
    dplyr::ungroup()

  # Calculate radii (Y-axis scaling)
  max_star_y <- max(pcoa_stats$SE_Max, pcoa_stats$Q75_Position, pcoa_stats$Mean_Position, na.rm = TRUE)
  if (is.infinite(max_star_y) || is.na(max_star_y)) max_star_y <- max(pcoa_stats$Mean_Position, na.rm = TRUE)
  
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

  # Convert geom_rect data into explicit polygons to avoid coord_radar() conversion issues
  iris_poly_list <- lapply(seq_len(nrow(iris_df)), function(i) {
    row <- iris_df[i, , drop = FALSE]
    data.frame(
      Sample = rep(row$Sample, 4),
      Taxa_Group = rep(row$Taxa_Group, 4),
      sample_var_col = rep(row[[sample_var]], 4),
      poly_x = c(row$xmin, row$xmin, row$xmax, row$xmax),
      poly_y = c(row$ymin, row$ymax, row$ymax, row$ymin),
      poly_group = rep(paste0(row$Sample, "_", row$Taxa_Group, "_", i), 4),
      stringsAsFactors = FALSE
    )
  })
  iris_poly_df <- do.call(rbind, iris_poly_list)
  colnames(iris_poly_df)[3] <- sample_var

  # --- 4. Plot Ordering ---
  if (!is.null(plot_order)) {
    pcoa_stats[[sample_var]] <- factor(pcoa_stats[[sample_var]], levels = plot_order)
    iris_poly_df[[sample_var]] <- factor(iris_poly_df[[sample_var]], levels = plot_order)
  }

  # --- 5. Colors ---
  groups <- unique(as.character(pcoa_stats[[sample_var]]))
  if (is.null(colors_star)) {
    if (exists("get_default_colors", mode = "function")) {
      colors_star <- get_default_colors(groups)
    } else {
      colors_star <- stats::setNames(grDevices::hcl.colors(length(groups), "Viridis"), groups)
    }
  }

  taxa_groups_all <- unique(iris_poly_df$Taxa_Group)
  if (is.null(colors_taxa)) {
    colors_taxa <- stats::setNames(grDevices::hcl.colors(length(taxa_groups_all), "Spectral"), taxa_groups_all)
  }

  # --- 6. Plot Construction ---
  title_suffix <- switch(error_bar, "IQR" = " with IQR", "SE" = " with SE", "none" = "")

  p <- ggplot2::ggplot() +
    # Layer 1: PCoA Star Polygon
    ggplot2::geom_polygon(data = pcoa_stats, ggplot2::aes(x = Axis_numeric, y = Mean_Position, group = !!dplyr::sym(sample_var), fill = !!dplyr::sym(sample_var), color = !!dplyr::sym(sample_var)), linewidth = 0.8, alpha = fill_alpha)
  
  # Error bars for PCoA Star
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
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%.2f", x - y_offset), breaks = scales::breaks_pretty(n = 3)) +
    coord_radar(clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "serif"), 
      axis.text.x = ggplot2::element_text(color = "black", size = 10),
      axis.text.y = ggplot2::element_text(color = "gray50"),
      panel.grid.major = ggplot2::element_line(color = "#e8e8e8", linewidth = 0.5),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    ) +
    ggplot2::labs(
      title = paste0("Mean ", method, " Position", title_suffix, " & Taxa Iris"),
      x = "",
      y = "Mean Score / Taxonomy"
    )

  if (view_type == "separate") {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", sample_var)))
  }

  return(p)
}
