#' Plot a Star Plot of Relative Abundance for Three Taxa
#'
#' This function generates a star (radar) plot showing the relative abundance
#' of three specified taxa plus an "Other" category, faceted by a sample variable.
#'
#' @param physeq A phyloseq object containing your microbiome data.
#' @param sample_var A character string. The name of the column in your sample data
#'   to use for faceting (e.g., "location", "treatment").
#' @param taxa_rank A character string. The taxonomic rank you are providing taxa
#'   names for (e.g., "Phylum", "Genus").
#' @param taxa_names A character vector containing the names of the three specific
#'   taxa you want to plot. If NULL (default), the top 4 most abundant taxa are selected.
#'
#' @param error_bar A character string, one of "IQR" (default), "SE", or "none".
#'   "IQR" draws the interquartile range of the samples in the group, around the
#'   group **median**; "SE" draws one standard error of the group **mean**,
#'   around the mean. The two describe different things — spread of the samples
#'   versus precision of the centre — and the plot title records which was used.
#' @param samplecolumn This is the ID column for your samples.
#' @param view_type A character string, either "separate" (default) or "together".
#'   Determines whether groups are plotted on the same plot or faceted.
#' @param fill_alpha A numeric value between 0 and 1 (default 0.2).
#'   Controls the transparency of the polygon fill under the star plot.
#' @param log_scale A logical value. If TRUE, applies a pseudo-log transformation to the y-axis.
#' @param plot_order A character vector for custom ordering of the sample variable,
#'   "hclust" for complete-linkage clustering on Bray-Curtis distances between the
#'   group mean abundances, or NULL (default) for alphabetical.
#' @param base_colors Optional. A character vector of base colors to use for hclust coloring.
#' @param colors_all Optional. A character vector of colors, or "hclust" for automatic clustering colors.
#' @param show_caption Logical. If FALSE (default), omits the explanatory caption below the plot.
#' @param axis_label_pos Character. Position of concentric-circle tick labels:
#'   `"side"` (default; along the left margin), `"middle"` (at the 12 o'clock
#'   position on each concentric circle), or `"none"` (hidden).
#' @param font Optional character string specifying the font family for plot
#'   text (e.g., `"sans"`, `"serif"`). If `NULL` (default), uses the active
#'   ggplot2 theme default without overriding font family.
#' @param label_color Character string specifying the color for all plot text
#'   and labels (default `"black"`).
#'
#' @return A ggplot object representing the star plot.
#'
#' @section Interpreting the plot:
#' The abundances are compositional: relative abundances within a sample sum to
#' 1, so the spokes are not independent — a rise in one taxon mechanically
#' depresses the others, and the pooled "Other" wedge is usually the largest
#' single value. Group values are the mean (or median) of per-sample relative
#' abundances, which is not the same as the relative abundance of pooled counts.
#' The area enclosed by a polygon depends on the order of the spokes and is not
#' an interpretable quantity; spokes are ordered by decreasing overall abundance
#' with "Other" last.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @import vegan
#'
#' @importFrom stats median quantile sd
#'
#' @export
#'
#' @examples
#' data("cassiopea_ps")
#' plot_taxa_star(
#'   physeq = cassiopea_ps,
#'   sample_var = "Body_cav",
#'   taxa_rank = "Phylum",
#'   samplecolumn = "Sample_NAME",
#'   view_type = "separate",
#'   error_bar = "SE"
#' )
plot_taxa_star <- function(physeq, sample_var, taxa_rank = "OTU", taxa_names = NULL, colors_all, samplecolumn, view_type = "separate", error_bar = "IQR", fill_alpha = 0.2, log_scale = FALSE, plot_order = NULL, base_colors = NULL, show_caption = FALSE, axis_label_pos = c("side", "middle", "none"), font = NULL, label_color = "black") {
  # --- 1. Input Validation and Conversion ---

  # Check if input is a data frame and convert to phyloseq if necessary
  if (inherits(physeq, "data.frame")) {
    message("Input is a data frame.I will assume all nonnumeric columns are metadata and
            all numeric columns except sample_var are data.To avoid this,
            you can always build a ps object yourself. Converting to phyloseq object...")

    # Identify numeric columns as taxa, others as sample data
    # We assume taxa columns are numeric and sample data are not (or are specific known columns)
    # We will assume that strictly numeric columns that are NOT `sample_var` or `samplecolumn` are taxa.

    # Identify potential taxa columns (numeric)
    numeric_cols <- sapply(physeq, is.numeric)

    # Exclude known sample variables from taxa candidates
    known_sample_vars <- c(sample_var, samplecolumn)
    taxa_candidates <- names(physeq)[numeric_cols]
    taxa_cols <- setdiff(taxa_candidates, known_sample_vars)

    if (length(taxa_cols) == 0) {
      stop("Error: No numeric taxa columns found in the data frame after excluding sample variables.")
    }

    # Extract Taxa and Sample Data
    taxa_mat <- as.matrix(physeq[, taxa_cols, drop = FALSE])
    sample_dat <- physeq[, !names(physeq) %in% taxa_cols, drop = FALSE]

    # Create Phyloseq components
    # OTU Table: Taxa are columns in input DF -> taxa_are_rows = FALSE
    otu <- phyloseq::otu_table(taxa_mat, taxa_are_rows = FALSE)

    # Sample Data
    sdata <- phyloseq::sample_data(sample_dat)

    # Taxonomy Table
    # We need a dummy taxonomy table that maps the column names (Taxa) to the requested rank.
    # `taxa_rank` argument will be used as the rank name in the taxonomy table.
    # Create a matrix with 1 column, named `taxa_rank`
    tax_mat <- matrix(colnames(taxa_mat), ncol = 1)
    colnames(tax_mat) <- taxa_rank
    rownames(tax_mat) <- colnames(taxa_mat)

    tax <- phyloseq::tax_table(tax_mat)

    # Combine into phyloseq object
    physeq <- phyloseq::phyloseq(otu, sdata, tax)
  }

  if (!inherits(physeq, "phyloseq")) {
    stop("Error: 'physeq' must be a VALID phyloseq object or a data frame.")
  }

  # Check for NAs in the OTU table and remove offending samples
  otu_mat <- as(phyloseq::otu_table(physeq), "matrix")
  if (any(is.na(otu_mat))) {
    warning("NAs found. Removing offending samples.")
    if (phyloseq::taxa_are_rows(physeq)) {
      valid_samples <- colnames(otu_mat)[!apply(is.na(otu_mat), 2, any)]
    } else {
      valid_samples <- rownames(otu_mat)[!apply(is.na(otu_mat), 1, any)]
    }
    physeq <- phyloseq::prune_samples(valid_samples, physeq)
  }

  if (!sample_var %in% phyloseq::sample_variables(physeq)) {
    stop("Error: '", sample_var, "' is not a valid sample variable in the phyloseq object/data frame.")
  }

  if (!is.null(taxa_names) && length(taxa_names) > 8) {
    stop("Error: 'taxa_names' must be fewer than 8 taxa (think about visual relevance!).")
  }
  # Check if rank exists.
  # Note: For converted DFs, we created the rank column, so this should pass.
  # For existing phyloseq objects, this remains a valid check.
  if (!taxa_rank %in% phyloseq::rank_names(physeq) && taxa_rank != "OTU") {
    stop("Error: '", taxa_rank, "' is not a valid taxonomic rank in the phyloseq object.")
  }

  view_type <- match.arg(view_type, c("together", "separate"))
  error_bar <- match.arg(error_bar, c("IQR", "SE", "none"))
  axis_label_pos <- match.arg(axis_label_pos)

  # --- 2. Data Transformation ---

  # to_relative() guards against empty samples, which would otherwise become a
  # row of NaN that propagates into every downstream mean.
  physeq_rel <- phyloseq::transform_sample_counts(physeq, to_relative)

  if (taxa_rank != "OTU") {
    physeq_rel <- agg_taxa(physeq_rel, rank = taxa_rank)
  } else {
    physeq_rel <- physeq_rel
    message("Taxa not aggregated. Using OTUs.")
  }

  # Melt the data into a long format for ggplot2
  df <- phyloseq::psmelt(physeq_rel)

  # --- 2b. Default Taxa Selection ---
  if (is.null(taxa_names)) {
    # Calculate total abundance for each taxon at the specified rank
    top_taxa <- df %>%
      group_by(!!rlang::sym(taxa_rank)) %>%
      summarise(TotalAbundance = sum(Abundance)) %>%
      arrange(desc(TotalAbundance)) %>%
      slice_head(n = 4) %>%
      pull(!!rlang::sym(taxa_rank)) %>%
      as.character()

    taxa_names <- top_taxa
    message("No taxa_names provided. Defaulting to top 4 most abundant: ", paste(taxa_names, collapse = ", "))
  }

  # --- 3. Taxa Grouping ---

  # Group the taxa: the three specified taxa and "Other"
  df_grouped <- df %>%
    # Use the taxa_rank column to mutate the new variable
    mutate(Taxa_Group = if_else(get(taxa_rank) %in% taxa_names, as.character(get(taxa_rank)), "Other")) %>%
    # Group by the specified sample variable and the new Taxa_Group
    group_by(!!rlang::sym(sample_var), !!rlang::sym(samplecolumn), Taxa_Group) %>%
    # Calculate the mean abundance for each group
    summarise(sum_Abundance = sum(Abundance), .groups = "drop")

  df_grouped_2 <- df_grouped %>%
    # Group by the specified sample variable and the new Taxa_Group
    group_by(!!rlang::sym(sample_var), Taxa_Group) %>%
    # Calculate the centre and spread for each group
    summarise(
      N_samples = dplyr::n(),
      mean_Abundance = mean(sum_Abundance),
      median_Abundance = stats::median(sum_Abundance),
      # IQR Stats
      Q25_Abundance = quantile(sum_Abundance, 0.25),
      Q75_Abundance = quantile(sum_Abundance, 0.75),
      # SE Stats
      SE_Abundance = sd(sum_Abundance) / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    mutate(
      SE_Min = mean_Abundance - SE_Abundance,
      SE_Max = mean_Abundance + SE_Abundance
    )

  # The IQR brackets the median, so plot the median when the IQR is shown; the
  # SE describes the mean, so plot the mean when the SE is shown. Microbiome
  # abundances are strongly right-skewed, so the two centres differ materially.
  center_col <- if (error_bar == "IQR") "median_Abundance" else "mean_Abundance"
  df_grouped_2$Center_Abundance <- df_grouped_2[[center_col]]
  center_label <- if (error_bar == "IQR") "Median" else "Mean"

  # Groups of one have no spread to show; say so rather than silently dropping
  # the interval.
  if (error_bar != "none") {
    thin <- unique(as.character(df_grouped_2[[sample_var]][df_grouped_2$N_samples < 2]))
    if (length(thin) > 0) {
      warning(
        "No interval can be drawn for group(s) with a single sample: ",
        paste(thin, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # The enclosed area of a star depends on the order of its spokes, so order
  # them deliberately: by decreasing overall abundance, with "Other" last.
  taxa_by_abundance <- df_grouped_2 %>%
    dplyr::filter(Taxa_Group != "Other") %>%
    dplyr::group_by(Taxa_Group) %>%
    dplyr::summarise(total = sum(Center_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::pull(Taxa_Group) %>%
    as.character()

  df_grouped_2$Taxa_Group <- order_taxa_group(df_grouped_2$Taxa_Group, levels_order = taxa_by_abundance)

  # geom_polygon() joins vertices in ROW order, not in x-axis order. Setting the
  # factor levels alone re-labels the spokes without re-tracing the shape, so the
  # polygon crosses over itself. The rows must be sorted to match the levels.
  df_grouped_2 <- df_grouped_2 %>%
    dplyr::arrange(!!rlang::sym(sample_var), Taxa_Group)

  # --- 3b. Plot Ordering ---
  # Bray-Curtis is appropriate here: the values being clustered are non-negative
  # relative abundances. (It would not be appropriate on ordination scores.)
  cluster_groups <- function() {
    wide_df <- df_grouped_2 %>%
      dplyr::select(!!rlang::sym(sample_var), Taxa_Group, Center_Abundance) %>%
      tidyr::pivot_wider(names_from = Taxa_Group, values_from = Center_Abundance, values_fill = list(Center_Abundance = 0))

    dist_m <- vegan::vegdist(wide_df %>% dplyr::select(-!!rlang::sym(sample_var)), method = "bray")
    hc_res <- stats::hclust(dist_m, method = "complete")
    hc_res$labels <- as.character(wide_df[[sample_var]])
    hc_res
  }

  if (!is.null(plot_order)) {
    if (length(plot_order) == 1 && plot_order == "hclust") {
      message("Ordering groups by complete-linkage clustering on Bray-Curtis distances between group abundances...")
      hc <- cluster_groups()
      ordered_groups <- hc$labels[hc$order]

      df_grouped_2[[sample_var]] <- factor(df_grouped_2[[sample_var]], levels = ordered_groups)
    } else {
      # Custom order provided by user
      df_grouped_2[[sample_var]] <- factor(df_grouped_2[[sample_var]], levels = plot_order)
    }
  }

  # --- 4. Plotting ---

  # Handle default colors if missing
  if (missing(colors_all) || is.null(colors_all) || (length(colors_all) == 1 && colors_all == "hclust")) {
    if (is.factor(df_grouped_2[[sample_var]])) {
      groups <- levels(df_grouped_2[[sample_var]])
    } else {
      groups <- as.character(unique(df_grouped_2[[sample_var]]))
    }

    if (!missing(colors_all) && !is.null(colors_all) && colors_all == "hclust") {
      # Order colors by hierarchical clustering automatically
      hc_res <- cluster_groups()
      ordered_names <- hc_res$labels[hc_res$order]

      colors_all <- get_hclust_colors(hc_res, ordered_names, base_colors = base_colors)
    } else {
      colors_all <- get_default_colors(groups)
    }
  } else if (!is.null(plot_order) && length(plot_order) == 1 && plot_order == "hclust") {
    # If hclust + color options specified, treat colors as bases for hclust base_palette
    hc_res <- cluster_groups()
    ordered_names <- hc_res$labels[hc_res$order]

    colors_all <- get_hclust_colors(hc_res, ordered_names, base_colors = colors_all)
  }

  title_suffix <- switch(error_bar,
    "IQR" = " with IQR",
    "SE" = " with SE",
    "none" = ""
  )

  caption_text <- if (isTRUE(show_caption)) {
    switch(error_bar,
      "IQR" = "Bars span the interquartile range of samples within each group, around the group median.",
      "SE" = "Bars span one standard error of the group mean.",
      "none" = NULL
    )
  } else {
    NULL
  }

  # Create the star plot
  star_plot <- ggplot(df_grouped_2, aes(
    x = Taxa_Group, y = Center_Abundance, group = !!rlang::sym(sample_var),
    fill = !!rlang::sym(sample_var), color = !!rlang::sym(sample_var)
  )) +
    geom_polygon(aes(), linewidth = 0.8, show.legend = FALSE, alpha = fill_alpha) +
    scale_fill_manual(values = colors_all) +
    scale_color_manual(values = colors_all) +
    theme_minimal() +
    theme(
      text = element_text(color = label_color),
      axis.text = element_text(color = label_color),
      axis.text.x = element_text(color = label_color, size = 9),
      axis.text.y = element_text(color = label_color),
      axis.title = element_text(color = label_color),
      plot.title = element_text(color = label_color),
      strip.text = element_text(color = label_color),
      legend.text = element_text(color = label_color),
      legend.title = element_text(color = label_color),
      panel.grid.major = element_line(color = "#e8e8e8", linewidth = 0.5),
      panel.spacing = grid::unit(3.5, "lines"),
      plot.caption = element_text(color = label_color, hjust = 0, size = 8),
      plot.margin = margin(15, 30, 15, 30)
    ) +
    labs(
      title = paste0(center_label, " Relative Abundance of Key Taxa", title_suffix),
      x = "",
      y = paste0(center_label, " Relative Abundance"),
      color = sample_var,
      fill = sample_var,
      caption = caption_text
    ) +
    coord_radar(clip = "off")

  if (!is.null(font)) {
    star_plot <- star_plot + theme(text = element_text(family = font, color = label_color))
  }

  if (view_type == "separate") {
    star_plot <- star_plot + facet_wrap(as.formula(paste("~", sample_var)))
  }

  # Add Error Bars conditionally
  if (error_bar == "IQR") {
    star_plot <- star_plot +
      geom_errorbar(aes(ymin = Q25_Abundance, ymax = Q75_Abundance), width = 0.2, show.legend = FALSE)
  } else if (error_bar == "SE") {
    star_plot <- star_plot +
      geom_errorbar(aes(ymin = SE_Min, ymax = SE_Max), width = 0.2, show.legend = FALSE)
  }

  if (log_scale) {
    star_plot <- star_plot + scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.01),
      breaks = scales::breaks_pretty(n = 3),
      expand = expansion(mult = c(0, 0.12))
    )
  } else {
    star_plot <- star_plot + scale_y_continuous(
      breaks = scales::breaks_pretty(n = 3),
      expand = expansion(mult = c(0, 0.12))
    )
  }

  star_plot <- apply_radial_labels(
    star_plot,
    axis_label_pos = axis_label_pos,
    label_fn = function(x) format(x, trim = TRUE),
    font = font,
    label_color = label_color
  )

  return(star_plot)
}
