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
#'   Determines the type of error bars to display.
#' @param samplecolumn This is the ID column for your samples.
#' @param view_type A character string, either "separate" (default) or "together".
#'   Determines whether groups are plotted on the same plot or faceted.
#' @param fill_alpha A numeric value between 0 and 1 (default 0.2).
#'   Controls the transparency of the polygon fill under the star plot.
#' @param log_scale A logical value. If TRUE, applies a pseudo-log transformation to the y-axis.
#' @param plot_order A character vector for custom ordering of the sample variable, "hclust" for Ward's clustering based on Euclidean distance, or NULL (default) for alphabetical.
#'
#' @return A ggplot object representing the star plot.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @import vegan
#'
#' @export
#'
#' @examples
#' # ...
#' #   plot_taxa_star(
#' #     physeq = gp_subset,
#' #     sample_var = "SampleType",
#' #     taxa_rank = "Phylum",
#' #     # taxa_names = c("Firmicutes", "Bacteroidetes"), # Optional
#' #     samplecolumn = "Sample_ID",
#' #     view_type = "separate",
#' #     error_bar = "SE",
#' #     fill_alpha = 0.2,
#' #     log_scale = FALSE,
#' #     plot_order = NULL
#' #    )
#' #
#'
#' # Example with a Data Frame
#' # ...
#' # plot_taxa_star(
#' #   ...
#' #   error_bar = "none"
#' # )
plot_taxa_star <- function(physeq, sample_var, taxa_rank = "OTU", taxa_names = NULL, colors_all, samplecolumn, view_type = "separate", error_bar = "IQR", fill_alpha = 0.2, log_scale = FALSE, plot_order = NULL) {
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

  # --- 2. Data Transformation ---

  physeq_rel <- phyloseq::transform_sample_counts(physeq, function(x) 1 * x / sum(x))

  if (taxa_rank != "OTU") {
    physeq_taxglom <- tax_agg(physeq_rel, rank = taxa_rank)
    physeq_rel <- physeq_taxglom
    message("Taxa aggregated to rank: ", taxa_rank, "using microViz::tax_agg.")
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
      group_by(!!sym(taxa_rank)) %>%
      summarise(TotalAbundance = sum(Abundance)) %>%
      arrange(desc(TotalAbundance)) %>%
      slice_head(n = 4) %>%
      pull(!!sym(taxa_rank)) %>%
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
    group_by(!!sym(sample_var), !!sym(samplecolumn), Taxa_Group) %>%
    # Calculate the mean abundance for each group
    summarise(sum_Abundance = sum(Abundance), .groups = "drop")

  df_grouped_2 <- df_grouped %>%
    # Group by the specified sample variable and the new Taxa_Group
    group_by(!!sym(sample_var), Taxa_Group) %>%
    # Calculate the mean abundance and stats for each group
    summarise(
      mean_Abundance = mean(sum_Abundance),
      # IQR Stats
      Q25_Abundance = quantile(sum_Abundance, 0.25),
      Q75_Abundance = quantile(sum_Abundance, 0.75),
      # SE Stats
      SE_Abundance = sd(sum_Abundance) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      SE_Min = mean_Abundance - SE_Abundance,
      SE_Max = mean_Abundance + SE_Abundance
    )

  # --- 3b. Plot Ordering ---
  if (!is.null(plot_order)) {
    if (length(plot_order) == 1 && plot_order == "hclust") {
      message("Calculating hclust for plot ordering using Bray-Curtis distance and complete linkage...")
      # Calculate hclust
      wide_df <- df_grouped_2 %>%
        dplyr::select(!!sym(sample_var), Taxa_Group, mean_Abundance) %>%
        tidyr::pivot_wider(names_from = Taxa_Group, values_from = mean_Abundance, values_fill = list(mean_Abundance = 0))

      dist_mat <- vegan::vegdist(wide_df %>% dplyr::select(-!!sym(sample_var)), method = "bray")
      hc <- stats::hclust(dist_mat, method = "complete")
      ordered_groups <- wide_df[[sample_var]][hc$order]

      df_grouped_2[[sample_var]] <- factor(df_grouped_2[[sample_var]], levels = ordered_groups)
    } else {
      # Custom order provided by user
      df_grouped_2[[sample_var]] <- factor(df_grouped_2[[sample_var]], levels = plot_order)
    }
  }

  # --- 4. Plotting ---

  # Handle default colors if missing
  if (missing(colors_all) || (length(colors_all) == 1 && colors_all == "hclust")) {
    if (is.factor(df_grouped_2[[sample_var]])) {
      groups <- levels(df_grouped_2[[sample_var]])
    } else {
      groups <- as.character(unique(df_grouped_2[[sample_var]]))
    }

    if (!missing(colors_all) && colors_all == "hclust") {
      # Order colors by hierarchical clustering automatically
      long_df <- df_grouped_2 %>%
        dplyr::select(!!sym(sample_var), Taxa_Group, mean_Abundance) %>%
        tidyr::pivot_wider(names_from = Taxa_Group, values_from = mean_Abundance, values_fill = list(mean_Abundance = 0))

      dist_m <- vegan::vegdist(long_df %>% dplyr::select(-!!sym(sample_var)), method = "bray")
      hc_res <- stats::hclust(dist_m, method = "complete")
      hc_res$labels <- as.character(long_df[[sample_var]])
      ordered_names <- hc_res$labels[hc_res$order]

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

  # Create the star plot
  star_plot <- ggplot(df_grouped_2, aes(
    x = Taxa_Group, y = mean_Abundance, group = !!sym(sample_var),
    fill = !!sym(sample_var), color = !!sym(sample_var)
  )) +
    geom_polygon(aes(), linewidth = 0.8, show.legend = FALSE, alpha = fill_alpha) +
    scale_fill_manual(values = colors_all) +
    scale_color_manual(values = colors_all) +
    theme_minimal() +
    theme(
      text = element_text(family = "serif"), axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "gray50"),
      panel.grid.major = element_line(color = "#e8e8e8", linewidth = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(
      title = paste0("Relative Abundance of Key Taxa", title_suffix),
      x = "",
      y = "Relative Abundance",
      color = "Sample"
    ) +
    coord_radar(clip = "off")

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
    star_plot <- star_plot + scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.01), breaks = scales::breaks_pretty(n = 3))
  } else {
    star_plot <- star_plot + scale_y_continuous(breaks = scales::breaks_pretty(n = 3))
  }

  return(star_plot)
}
