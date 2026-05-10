# ordi_star_hclust_tree.R
# One-time script to generate the hclust tree underlying the plot_ordi_star 
# command, with colored dots on the tips that match the plot output.

library(phyloseq)
library(dplyr)
library(tidyr)
library(vegan)
library(stats)

# Source the colors utility to load get_hclust_colors if it isn't already loaded
if (!exists("get_hclust_colors", mode = "function")) {
  tryCatch({
    source("R/utils-colors.R")
  }, error = function(e) {
    warning("Could not find 'R/utils-colors.R'. Please make sure your working directory is set to the package root, or that 'get_hclust_colors' is already loaded.")
  })
}

#' Generate hclust tree with matching plot_ordi_star colors
#' 
#' @param physeq A phyloseq object
#' @param sample_var The sample variable used for grouping
#' @param method Ordination method (default "PCoA")
#' @param distance Distance metric (default "bray")
#' @param n_axes Number of axes to calculate the tree on (default 5, matches plot_ordi_star)
#' @param base_colors Optional vector of base colors to pass to get_hclust_colors
#' 
plot_ordi_star_tree <- function(physeq, sample_var, ord = NULL, method = "PCoA", distance = NULL, n_axes = 5, base_colors = NULL, colors_all = NULL) {
  
  # 1. Ordination (Use existing if provided)
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
  
  # Extract axes
  if (inherits(ord, c("pcoa", "dpcoa", "list")) && !is.null(ord$vectors)) {
    vecs <- ord$vectors
    axis_raw_names <- paste0("Axis ", 1:ncol(vecs))
  } else if (inherits(ord, c("metaMDS", "monoMDS"))) {
    vecs <- ord$points
    axis_raw_names <- paste0("NMDS", 1:ncol(vecs))
  } else if (inherits(ord, c("rda", "cca", "prcomp"))) {
    vecs <- vegan::scores(ord, display = "sites")
    axis_raw_names <- paste0("PC", 1:ncol(vecs))
  } else {
    vecs <- tryCatch(vegan::scores(ord, display = "sites"), error = function(e) ord$vectors)
    axis_raw_names <- paste0("Axis ", 1:ncol(vecs))
  }
  
  actual_axes <- min(n_axes, ncol(vecs))
  axis_labels <- axis_raw_names[1:actual_axes]
  
  pcoa_axes <- data.frame(vecs[, 1:actual_axes, drop = FALSE])
  colnames(pcoa_axes) <- axis_labels
  
  sample_dat <- data.frame(phyloseq::sample_data(physeq))
  pcoa_df <- cbind(pcoa_axes, sample_dat)
  
  # Apply identical shift
  min_val <- min(pcoa_axes, na.rm = TRUE)
  if (min_val < 0) {
    y_offset <- abs(min_val) * 1.1
  } else {
    y_offset <- 0
  }
  
  pcoa_df_shifted <- pcoa_df
  pcoa_df_shifted[axis_labels] <- pcoa_df[axis_labels] + y_offset
  
  pcoa_long_raw <- pcoa_df_shifted %>%
    tidyr::pivot_longer(cols = all_of(axis_labels), names_to = "Axis", values_to = "Value")
  
  pcoa_stats <- pcoa_long_raw %>%
    group_by(!!sym(sample_var), Axis) %>%
    summarise(
      Mean_Position = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2. Calculate hclust exactly as in plot_ordi_star
  wide_df <- pcoa_stats %>%
    dplyr::select(!!sym(sample_var), Axis, Mean_Position) %>%
    tidyr::pivot_wider(names_from = Axis, values_from = Mean_Position, values_fill = list(Mean_Position = 0))
  
  dist_mat <- vegan::vegdist(wide_df %>% dplyr::select(-!!sym(sample_var)), method = "bray")
  hc_res <- stats::hclust(dist_mat, method = "complete")
  
  # Crucial: Apply labels to the hclust object so the dendrogram has names
  hc_res$labels <- as.character(wide_df[[sample_var]])
  
  ordered_names <- wide_df[[sample_var]][hc_res$order]
  
  # 3. Get corresponding hclust colors
  if (is.null(colors_all) || (length(colors_all) == 1 && colors_all == "hclust")) {
    tree_colors <- get_hclust_colors(hc_res, ordered_names, base_colors = base_colors)
  } else {
    # If colors_all is provided, it acts as the base_colors just like in plot_ordi_star
    tree_colors <- get_hclust_colors(hc_res, ordered_names, base_colors = colors_all)
  }
  
  # 4. Create dendrogram and add colored dots
  dend <- as.dendrogram(hc_res)
  
  # Custom dendrapply function to set node parameters for leaves
  color_leaves <- function(node, colors) {
    if (is.leaf(node)) {
      leaf_name <- as.character(attr(node, "label"))
      
      # Safely extract color (avoids subscript out of bounds if leaf_name was numeric or missing)
      if (!is.null(names(colors)) && leaf_name %in% names(colors)) {
        leaf_color <- colors[[leaf_name]]
      } else if (is.null(names(colors)) && length(colors) > 0) {
        leaf_color <- colors[1] # Fallback for unnamed vectors
      } else {
        leaf_color <- "black"   # Fallback if name not found
      }
      
      # Set plotting character (pch = 19 is a solid circle) and color
      attr(node, "nodePar") <- c(attr(node, "nodePar"), list(pch = 19, col = leaf_color, cex = 2))
    }
    return(node)
  }
  
  dend_colored <- dendrapply(dend, color_leaves, colors = tree_colors)
  
  # 5. Plot
  plot(dend_colored, 
       main = paste("Hierarchical Clustering of", sample_var, "\n(Matches plot_ordi_star colors)"), 
       ylab = "Bray-Curtis Distance")
  
  # Return data invisibly
  invisible(list(hc = hc_res, dend = dend_colored, colors = tree_colors))
}

# ==============================================================================
# HOW TO USE THIS SCRIPT:
# 
# 1. Load your phyloseq object into your R environment (e.g., 'my_physeq').
# 2. Run the function below, replacing the arguments with your actual data.
# 
# Example:
# plot_ordi_star_tree(
#   physeq = my_physeq,
#   sample_var = "YourSampleGroupingColumn",
#   method = "PCoA", 
#   distance = "bray",
#   n_axes = 5
# )
# ==============================================================================
