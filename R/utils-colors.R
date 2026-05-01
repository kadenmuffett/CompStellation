#' Get Default Muted Colors
#'
#' Returns a named vector of colors based on the number of groups.
#' Uses a muted, pleasing palette (Okabe-Ito inspired or similar).
#'
#' @param groups A character vector of group names.
#'
#' @return A named character vector associated with the groups.
#' @keywords internal
get_default_colors <- function(groups) {
  n <- length(groups)

  # A pleasing muted palette (Okabe-Ito inspired + extras)
  # 1. Orange, 2. Sky Blue, 3. Bluish Green, 4. Yellow, 5. Blue, 6. Vermilion, 7. Reddish Purple, 8. Grey
  base_palette <- c(
    "#99cf4d", "#e7c03f", "#2206b0", "#4b070f", "#d115c7", "#d474e0", "#9204c6", "#484ab1",
    "#bb206e", "#AA4499", "#332288", "#e7d399", "#117733", "#44AA99"
  )

  if (n > length(base_palette)) {
    warning("More groups than default colors available. Recycling colors.")
    colors <- rep(base_palette, length.out = n)
  } else {
    colors <- base_palette[1:n]
  }

  names(colors) <- groups
  return(colors)
}

get_default_colors_hclust <- function(groups) {
  n <- length(groups)

  # A pleasing muted palette (Okabe-Ito inspired + extras)
  # 1. Orange, 2. Sky Blue, 3. Bluish Green, 4. Yellow, 5. Blue, 6. Vermilion, 7. Reddish Purple, 8. Grey
  base_palette <- c("#990000", # Crimson Red
  "#994C00", # Burnt Orange
  "#666600", # Dark Citrine
  "#006600", # Deep Kelly Green
  "#006666", # Dark Manganese Teal
  "#000099", # Royal Blue
  "#4C0099", # Deep Violet
  "#99004C"  # Raspberry Red)


  if (n > length(base_palette)) {
    warning("More groups than default colors available. Recycling colors.")
    colors <- rep(base_palette, length.out = n)
  } else {
    colors <- base_palette[1:n]
  }

  names(colors) <- groups
  return(colors)
}

#' Get Clade-Based Colors from Hierarchical Clustering
#'
#' Cuts a dendrogram into clades and assigns a visually distinct base color
#' to each clade. Leaves within the same clade are assigned shades of that base color.
#'
#' @param hc_res An hclust object.
#' @param ordered_names The names of the elements in the order of the clustered leaves.
#' @param k Optional. The number of clades to cut the tree into. If NULL, determines dynamically.
#'
#' @return A named character vector of hex colors matching ordered_names.
#' @keywords internal
get_hclust_colors <- function(hc_res, ordered_names, k = NULL) {
  ordered_names <- as.character(ordered_names)
  n <- length(ordered_names)
  if (n <= 3) {
    return(stats::setNames(get_default_colors_hclust(as.character(1:n)), ordered_names))
  }

  if (is.null(k)) {
    # Dynamically pick a sensible number of clades
    k <- max(2, min(5, ceiling(n / 3)))
  }
  k <- min(k, n) # prevent k > n

  # Reconstruct original labels to ensure cutree uses proper names
  original_names <- as.character(ordered_names[order(hc_res$order)])
  hc_res$labels <- original_names

  clusters <- stats::cutree(hc_res, k = k)
  base_colors <- get_default_colors_hclust(as.character(1:k))

  final_colors <- character(n)
  names(final_colors) <- ordered_names

  for (i in 1:k) {
    clade_members <- as.character(names(clusters[clusters == i]))
    ordered_clade_members <- intersect(ordered_names, clade_members)
    n_members <- length(ordered_clade_members)

    if (n_members == 1) {
      final_colors[ordered_clade_members] <- as.character(base_colors[i])
    } else if (n_members > 1) {
      # Mix the base color with white to get distinct shades within the clade
      ramp <- grDevices::colorRampPalette(c("#FFFFFF", base_colors[i]))
      shades <- ramp(n_members + 2)[3:(n_members + 2)]

      # Assign shades according to the hclust tree leaf order
      final_colors[ordered_clade_members] <- shades
    }
  }
  # Fallback to prevent ggplot exceptions
  final_colors[final_colors == "" | is.na(final_colors)] <- "#999999"

  return(final_colors)
}
