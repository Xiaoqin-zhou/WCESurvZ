#' plot.WCEKM
#'
#' Plot Weighted Composite Endpoint Kaplan-Meier Survival Curves
#'
#' This function generates a Kaplan-Meier survival plot with weighted composite endpoints from a WCEKM object.
#'
#' @param WCEKM_obj An object of class "WCEKM" containing the Kaplan-Meier survival data.
#' @param xlim The limits for the x-axis. Default is NULL, which determines the range automatically.
#' @param ylim The limits for the y-axis. Default is c(0, 1).
#' @param xlab The label for the x-axis. Default is NULL.
#' @param ylab The label for the y-axis. Default is NULL.
#' @param x_ticks The major ticks for the x-axis. Default is NULL.
#' @param y_ticks The major ticks for the y-axis. Default is NULL.
#' @param base_size The base size for the theme. Default is 12.
#' @param title The title of the plot. Default is NULL.
#' @param colors A vector of colors for the different groups. Default is NULL.
#' @param legend_position The position of the legend. Default is c(0.9, 0.8).
#' @param legend_title The title of the legend. Default is NULL.
#' @param legend_labels The labels for the legend. Default is NULL.
#' @param conf.int.alpha The alpha transparency for the confidence intervals. Default is 0.3.
#' @param theme The ggplot2 theme to use. Default is theme_WCESurvZ(base_size = base_size).
#' @param ... Additional arguments passed to methods.
#'
#' @return This function does not return a value. It generates a plot.
#' @export
#'
#' @examples
#' # Example usage:
#' weights <- c(CVdeath = 1, MI = 0.55, Stroke = 0.455)
#' WCE_obj <- WCE_KMSurv(surv_exdata, weights)
#' plot(WCE_obj)
#'
plot.WCEKM <- function(WCEKM_obj, xlim = NULL, ylim = c(0, 1),
                       xlab = NULL, ylab = NULL,  x_ticks = NULL, y_ticks = NULL,
                       base_size = 12,title = NULL, pval.coord = c(Inf,Inf), pval.size = 6,
                       colors = NULL, legend_position = c(0.9, 0.8),
                       legend_title = NULL, legend_labels = NULL,conf.int.alpha = 0.3,theme = theme_WCESurvZ(base_size = base_size), ...) {

  # Extract km_data from the WCEKM object
  km_data <- WCEKM_obj$km_data

  # If x_limits is NULL, determine the x-axis range automatically
  if (is.null(xlim)) {
    xlim <- c(0, max(km_data$time, na.rm = TRUE))
  }

  # If y_ticks is NULL, set the y-axis ticks to 0.1 increments
  if (is.null(y_ticks)) {
    y_ticks <- seq(ylim[1], ylim[2], by = 0.1)
  }

  # If x_ticks is NULL, determine the x-axis ticks automatically
  if (is.null(x_ticks)) {
    max_time <- xlim[2]
    digit_count <- 10^(floor(log10(max_time)))
    rounded_tick <- ceiling(max_time / digit_count)
    x_ticks <- seq(xlim[1], xlim[2], by = rounded_tick * digit_count / 10)
  }

  # Set default x-axis and y-axis labels
  if (is.null(xlab)) {
    xlab <- "time"
  }
  if (is.null(ylab)) {
    ylab <- "Survival probability"
  }

  # Set default legend title
  if (is.null(legend_title)) {
    legend_title <- "strata"
  }

  # Set default legend labels
  if (is.null(legend_labels)) {
    legend_labels <- unique(km_data$group)
  }

  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(WCEKM_obj$groups))
  }

  # Use ggplot2 to plot weighted Kaplan-Meier survival curves
  library(ggplot2)


  p_value_formatted <- if (WCEKM_obj$pvalue < 0.0001) {
    formatC(WCEKM_obj$pvalue, format = "e", digits = 2)
  } else {
    sprintf("%.4f", WCEKM_obj$pvalue)
  }

  # Create the Kaplan-Meier plot
  p <- ggplot(km_data, aes(x = time, y = survival,color = group, group = group)) +
    geom_step(size = 1.02) +
    labs(title = title, x = xlab, y = ylab) +
    geom_point(data = subset(km_data, n.censor > 0), aes(x = time, y = survival), shape = 3, size = 3) +
    theme +
    geom_confint(mapping = aes(ymin = lower, ymax = upper, fill = group), data = km_data, stat = "identity",
        position = "identity", na.rm = TRUE,alpha = conf.int.alpha, color = NA)+
    theme(legend.position = legend_position) +
    annotate("text", x = pval.coord[1], y = pval.coord[2], label = paste("p =", p_value_formatted), hjust = 1.5, vjust = 2, size = pval.size) +
    scale_color_manual(values = colors, labels = legend_labels) +
    scale_fill_manual(values = colors, labels = legend_labels) +
    scale_x_continuous(
      breaks = x_ticks,  # Set major ticks
      labels = as.character(x_ticks),  # Custom tick labels
      limits = xlim,  # Set x-axis limits
      expand = c(0.01, 0),  # No extra space
      position = "bottom"  # Place x-axis at the bottom
    ) +
    scale_y_continuous(
      breaks = y_ticks,  # Set major ticks
      limits = ylim  # Set y-axis limits
    )


  print(p)
}
