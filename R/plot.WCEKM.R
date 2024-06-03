#' Plot Weighted Composite Endpoint Kaplan-Meier Survival Curves
#'
#' This function plots Kaplan-Meier survival curves for weighted composite endpoints.
#'
#' @param WCEKM_obj An object of class "WCEKM" containing Kaplan-Meier survival data.
#' @param x_limits A vector specifying the limits for the x-axis. Default is NULL.
#' @param y_limits A vector specifying the limits for the y-axis. Default is c(0, 1).
#' @param font_size The base font size for the plot. Default is 12.
#' @param title The title of the plot. Default is "非加权Kaplan-Meier生存曲线".
#' @param x_label The label for the x-axis. Default is NULL.
#' @param y_label The label for the y-axis. Default is NULL.
#' @param legend_position A vector specifying the position of the legend. Default is c(0.9, 0.8).
#' @param colors A vector of colors for the different groups. Default is NULL.
#' @param x_ticks A vector specifying the ticks for the x-axis. Default is NULL.
#' @param y_ticks A vector specifying the ticks for the y-axis. Default is NULL.
#' @param legend_title The title for the legend. Default is NULL.
#' @param legend_labels A vector of labels for the legend. Default is NULL.
#'
#' @return A ggplot2 object showing the Kaplan-Meier survival curves.
#' @export
#'
#' @examples
#' # Example usage:
#' weights <- c(CVdeath = 1, MI = 0.55, Stroke = 0.455)
#' WCE_obj <- WCE_KMSurv(surv_exdata, weights)
#' plot(WCE_obj)
#'
plot.WCEKM <- function(WCEKM_obj, x_limits = NULL, y_limits = c(0, 1),
                       title = NULL,base_size = 12,
                       x_label = NULL, y_label = NULL,  x_ticks = NULL, y_ticks = NULL,
                       colors = NULL, legend_position = c(0.9, 0.8),
                       legend_title = NULL, legend_labels = NULL,theme = theme_WCESurvZ(base_size = base_size)) {

  # Extract km_data from the WCEKM object
  km_data <- WCEKM_obj$km_data

  # If x_limits is NULL, determine the x-axis range automatically
  if (is.null(x_limits)) {
    x_limits <- c(0, max(km_data$time, na.rm = TRUE))
  }

  # If y_ticks is NULL, set the y-axis ticks to 0.1 increments
  if (is.null(y_ticks)) {
    y_ticks <- seq(y_limits[1], y_limits[2], by = 0.1)
  }

  # If x_ticks is NULL, determine the x-axis ticks automatically
  if (is.null(x_ticks)) {
    max_time <- x_limits[2]
    digit_count <- 10^(floor(log10(max_time)))
    rounded_tick <- ceiling(max_time / digit_count)
    x_ticks <- seq(x_limits[1], x_limits[2], by = rounded_tick * digit_count / 10)
  }

  # Set default x-axis and y-axis labels
  if (is.null(x_label)) {
    x_label <- "time"
  }
  if (is.null(y_label)) {
    y_label <- "Survival probability"
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
    geom_step() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, color = NA) +
    labs(title = title, x = x_label, y = y_label) +
    geom_point(data = subset(km_data, n.censor > 0), aes(x = time, y = survival), shape = 3, size = 2) +
    theme +
    theme(legend.position = legend_position) +
    annotate("text", x = Inf, y = Inf, label = paste("P-value:", p_value_formatted), hjust = 1.5, vjust = 2, size = base_size/3) +
    scale_color_manual(values = colors, labels = legend_labels) +
    scale_fill_manual(values = colors, labels = legend_labels) +
    scale_x_continuous(
      breaks = x_ticks,  # Set major ticks
      labels = as.character(x_ticks),  # Custom tick labels
      limits = x_limits,  # Set x-axis limits
      expand = c(0.01, 0),  # No extra space
      position = "bottom"  # Place x-axis at the bottom
    ) +
    scale_y_continuous(
      breaks = y_ticks,  # Set major ticks
      limits = y_limits  # Set y-axis limits
    )


  # Print the plot
  print(p)
}
