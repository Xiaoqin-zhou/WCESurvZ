#' Connect observations by stairs.
#'
#' @param mapping the aesthetic mapping
#' @param data a layer specific dataset
#' @param stat the statistical transformation to use on the data for this layer
#' @param position the position adjustment to use for overlapping points on this layer
#' @param na.rm logical frag whether silently remove missing values
#' @param ... other arguments passed to methods
geom_confint <- function (mapping = NULL, data = NULL, stat = "identity",
                          position = "identity", na.rm = FALSE, ...) {
  ggplot2::layer(mapping = mapping,
                 data = data,
                 stat = stat,
                 geom = GeomConfint,
                 position = position,
                 params = list(na.rm = na.rm, ...))
}

GeomConfint <- ggplot2::ggproto('GeomConfint', ggplot2::GeomRibbon,
  required_aes = c("x", "ymin", "ymax"),

  # Function to draw the confidence intervals for a group
  draw_group = function(self, data, panel_scales, coord, na.rm = FALSE) {
    # If na.rm is TRUE, remove rows with missing required aesthetic values
    if (na.rm) data <- data[stats::complete.cases(self$required_aes), ]

    # Order data by group and x
    data <- data[order(data$group, data$x), ]

    # Convert data into stairstep format for confidence intervals
    data <- self$stairstep_confint(data)

    # Call the draw_group method of the GeomRibbon class to draw the confidence intervals
    ggplot2:::GeomRibbon$draw_group(data, panel_scales, coord, na.rm = FALSE)
  },

  # Function to convert the data into a stairstep format for confidence intervals
  stairstep_confint = function (data) {
    # Convert data to data frame and order by x
    data <- as.data.frame(data)[order(data$x), ]

    # Number of rows in the data
    n <- nrow(data)

    # Create indices for y values in a stairstep format
    ys <- rep(1:n, each = 2)[-2 * n]

    # Create indices for x values in a stairstep format
    xs <- c(1, rep(2:n, each = 2))

    # Create a new data frame in stairstep format for confidence intervals
    data.frame(
      x = data$x[xs],
      ymin = data$ymin[ys],
      ymax = data$ymax[ys],
      data[xs, setdiff(names(data), c("x", "ymin", "ymax"))]
    )
  }
)
