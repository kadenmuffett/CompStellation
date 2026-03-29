#' Custom Radar Coordinate System
#' 
#' @param theta Variable to map angle to ('x' or 'y').
#' @param start Offset of starting point from 12 o'clock in radians.
#' @param direction 1, clockwise; -1, anticlockwise.
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
