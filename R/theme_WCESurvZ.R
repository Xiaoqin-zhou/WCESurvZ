#' theme_WCESurvZ
#'
#' Custom Theme for Weighted Composite Endpoint Survival Plots
#'
#' This function generates a custom ggplot2 theme tailored for weighted composite endpoint survival plots, with specific font settings for various plot elements.
#'
#' @param base_size Base font size for the theme. Default is 12.
#' @param base_family Base font family for the theme. Default is "".
#' @param font.main Font settings for the main title. Default is c(16, "plain", "black").
#' @param font.submain Font settings for the subtitle. Default is c(15, "plain", "black").
#' @param font.x Font settings for the x-axis title. Default is c(14, "plain", "black").
#' @param font.y Font settings for the y-axis title. Default is c(14, "plain", "black").
#' @param font.caption Font settings for the plot caption. Default is c(15, "plain", "black").
#' @param font.tickslab Font settings for the axis tick labels. Default is c(12, "plain", "black").
#' @param font.legend Font settings for the legend text. Default is c(10, "plain", "black").
#' @param ... Additional arguments passed to methods.
#'
#' @return A ggplot2 theme object customized for weighted composite endpoint survival plots.
#' @export
#'
#' @examples
#' # Example usage:
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_WCESurvZ()
#' print(p)
theme_WCESurvZ <- function(base_size = 12, base_family = "", font.main = c(16, "plain", "black"),
                           font.submain = c(15, "plain", "black"), font.x = c(14, "plain", "black"),
                           font.y = c(14, "plain", "black"), font.caption = c(15, "plain", "black"),
                           font.tickslab = c(12, "plain", "black"), font.legend = c(10, "plain", "black"),
                           ...) {

  # Parse font settings for various plot elements
  font.main <- .parse_font(font.main)
  font.x <- .parse_font(font.x)
  font.y <- .parse_font(font.y)
  font.submain <- .parse_font(font.submain)
  font.caption <- .parse_font(font.caption)
  font.tickslab <- .parse_font(font.tickslab)
  font.legend <- .parse_font(font.legend)

  # Define text elements for axis ticks and legend text
  tickslab <- element_text(size = font.tickslab$size, face = font.tickslab$face,
                           colour = font.tickslab$color, angle = 0)
  legend.text <- element_text(size = font.legend$size, face = font.legend$face,
                              colour = font.legend$color)

  # Create the custom theme based on theme_classic
  result <- theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = font.main$size, lineheight = 1, face = font.main$face,
                                colour = font.main$color),
      plot.subtitle = element_text(size = font.submain$size, lineheight = 1, face = font.submain$face,
                                   colour = font.submain$color),
      axis.title.x = element_text(size = font.x$size, face = font.x$face, colour = font.x$color),
      axis.title.y = element_text(angle = 90, size = font.y$size, face = font.y$face, colour = font.y$color),
      plot.caption = element_text(size = font.caption$size, lineheight = 1, face = font.caption$face,
                                  colour = font.caption$color),
      axis.text.x = tickslab,
      axis.text.y = tickslab,
      legend.text = legend.text,
      legend.title = legend.text
    )

  # Assign the "theme" class to the result
  class(result) <- "theme"

  return(result)
}




.parse_font <- function(font) {
  # Check if font is NULL
  if (is.null(font)) {
    res <- NULL
  }
  # Check if font inherits from list
  else if (inherits(font, "list")) {
    res <- font
  }
  # Parse font settings from character vector
  else {
    # Extract size using regular expression
    size <- grep("^[0-9]+$", font, perl = TRUE)
    face <- grep("plain|bold|italic|bold.italic", font, perl = TRUE)
    if (length(size) == 0)
      size <- NULL
    else
      size <- as.numeric(font[size])

    if (length(face) == 0)
      face <- NULL
    else
      face <- font[face]

    # Extract color settings
    color <- setdiff(font, c(size, face))
    if (length(color) == 0)
      color <- NULL

    # Create a list with parsed font settings
    res <- list(size = size, face = face, color = color)
  }
  res
}
