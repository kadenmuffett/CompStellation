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
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#882255", "#AA4499", "#332288", "#DDCC77", "#117733", "#44AA99"
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
  n <- length(ordered_names)
  if (n <= 3) {
    return(stats::setNames(get_default_colors(as.character(1:n)), ordered_names))
  }
  
  if (is.null(k)) {
    # Dynamically pick a sensible number of clades
    k <- max(2, min(5, ceiling(n / 2.5)))
  }
  k <- min(k, n) # prevent k > n
  
  clusters <- stats::cutree(hc_res, k = k)
  base_colors <- get_default_colors(as.character(1:k))
  
  final_colors <- character(n)
  names(final_colors) <- ordered_names
  
  for (i in 1:k) {
    clade_members <- names(clusters[clusters == i])
    n_members <- length(clade_members)
    
    if (n_members == 1) {
      final_colors[clade_members] <- base_colors[i]
    } else {
      # Mix the base color with white to get distinct shades within the clade
      ramp <- grDevices::colorRampPalette(c("#FFFFFF", base_colors[i]))
      shades <- ramp(n_members + 2)[3:(n_members + 2)]
      
      # Assign shades according to the hclust tree leaf order
      ordered_clade_members <- intersect(ordered_names, clade_members)
      final_colors[ordered_clade_members] <- shades
    }
  }
  
  return(final_colors)
}
