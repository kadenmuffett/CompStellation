#' Custom Radar Coordinate System
#' 
#' @param theta Variable to map angle to ('x' or 'y').
#' @param start Offset of starting point from 12 o'clock in radians.
#' @param direction 1, clockwise; -1, anticlockwise.
#' @param clip A character string for clipping (default "off").
#' 
#' @keywords internal
coord_radar <- function (theta = "x", start = 0, direction = 1, clip = "off")
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction), clip = clip,
          is_linear = function(coord) TRUE)
}

#' Check if coordinate system is linear (for radar plots)
#' 
#' @param coord Coordinate system object.
#' 
#' @keywords internal
is.linear.radar <- function(coord) TRUE

#' Apply radial axis label positioning to a radar ggplot object
#'
#' @param p A ggplot object using `coord_radar()`.
#' @param axis_label_pos One of `"side"` (default; standard left-side y-axis
#'   labels), `"middle"` (labels placed at 12 o'clock on each concentric ring),
#'   or `"none"` (radial tick labels hidden).
#' @param label_fn A function formatting numeric break values into character
#'   labels.
#' @param font Optional font family string. If `NULL`, inherits default font.
#' @param label_color Character string specifying the text color (default `"black"`).
#' @param r_max Optional numeric upper bound on radial breaks to label (used by
#'   `plot_ordi_star_iris()` so the outer iris ring is not given ordination
#'   score labels).
#'
#' @return The modified ggplot object.
#' @keywords internal
#' @noRd
apply_radial_labels <- function(p, axis_label_pos = c("side", "middle", "none"),
                                label_fn = function(x) format(x, trim = TRUE),
                                font = NULL, label_color = "black", r_max = Inf) {
  axis_label_pos <- match.arg(axis_label_pos)

  if (axis_label_pos == "side") {
    return(p)
  }

  if (axis_label_pos == "none") {
    return(p + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ))
  }

  # axis_label_pos == "middle": place labels at 12 o'clock (x = 0.5) on rings
  pb <- ggplot2::ggplot_build(p)
  y_scale <- pb$layout$panel_scales_y[[1]]
  breaks <- y_scale$get_breaks()
  breaks <- breaks[!is.na(breaks) & is.finite(breaks) & breaks <= r_max]

  if (length(breaks) > 0) {
    labels <- label_fn(breaks)
    keep <- !is.na(labels) & nzchar(trimws(labels))
    breaks <- breaks[keep]
    labels <- labels[keep]
  }

  p <- p + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

  if (length(breaks) > 0) {
    df_rings <- data.frame(.ring_y = breaks, .ring_label = labels, stringsAsFactors = FALSE)
    text_args <- list(
      data = df_rings,
      mapping = ggplot2::aes(x = 0.5, y = .ring_y, label = .ring_label),
      inherit.aes = FALSE,
      size = 2.7,
      color = label_color,
      vjust = -0.25,
      hjust = 0.5
    )
    if (!is.null(font)) {
      text_args$family <- font
    }
    p <- p + do.call(ggplot2::geom_text, text_args)
  }

  p
}
