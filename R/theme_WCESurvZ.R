#' theme_WCESurvZ
#'
#' @param base_size NA
#' @param base_family NA
#' @param font.main NA
#' @param font.submain NA
#' @param font.x NA
#' @param font.y NA
#' @param font.caption NA
#' @param font.tickslab NA
#' @param legend NA
#' @param font.legend NA
#' @param ... NA
#'
#' @return NA
#' @export
#'
#' @examples
theme_WCESurvZ <- function (base_size = 12, base_family = "", font.main = c(16,
    "plain", "black"), font.submain = c(15, "plain", "black"),
    font.x = c(14, "plain", "black"), font.y = c(14, "plain",
        "black"), font.caption = c(15, "plain", "black"), font.tickslab = c(12,
        "plain", "black"),  font.legend = c(10, "plain", "black"),
    ...)
{
    font.main <- .parse_font(font.main)
    font.x <- .parse_font(font.x)
    font.y <- .parse_font(font.y)
    font.submain <- .parse_font(font.submain)
    font.caption <- .parse_font(font.caption)
    font.tickslab <- .parse_font(font.tickslab)
    font.legend <- .parse_font(font.legend)

    tickslab <- element_text(size = font.tickslab$size, face = font.tickslab$face,
        colour = font.tickslab$color, angle = 0)
    legend.text <- element_text(size = font.legend$size, face = font.legend$face,
        colour = font.legend$color)

    result <- theme_classic(base_size = base_size, base_family = base_family) +
        theme(plot.title = element_text(size = font.main$size,
              lineheight = 1, face = font.main$face, colour = font.main$color),
              plot.subtitle = element_text(size = font.submain$size,
              lineheight = 1, face = font.submain$face, colour = font.submain$color),
              axis.title.x = element_text(size = font.x$size, face = font.x$face,
                colour = font.x$color),
              axis.title.y = element_text(angle = 90,
                size = font.y$size, face = font.y$face, colour = font.y$color),
              plot.caption = element_text(size = font.caption$size,
                lineheight = 1, face = font.caption$face, colour = font.caption$color),
              axis.text.x = tickslab,
              axis.text.y = tickslab,
              legend.text = legend.text, legend.title = legend.text)
    class(result) <- "theme"
    result
}



.parse_font <- function (font)
{
    if (is.null(font))
        res <- NULL
    else if (inherits(font, "list"))
        res <- font
    else {
        size <- grep("^[0-9]+$", font, perl = TRUE)
        face <- grep("plain|bold|italic|bold.italic", font, perl = TRUE)
        if (length(size) == 0)
            size <- NULL
        else size <- as.numeric(font[size])
        if (length(face) == 0)
            face <- NULL
        else face <- font[face]
        color <- setdiff(font, c(size, face))
        if (length(color) == 0)
            color <- NULL
        res <- list(size = size, face = face, color = color)
    }
    res
}
