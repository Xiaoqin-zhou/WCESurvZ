#' WCE_KMSurv
#'
#' Generate Weighted Composite Endpoint Kaplan-Meier Survival Data
#'
#' This function calculates weighted Kaplan-Meier survival curves for recurrent events in different groups.
#'
#' @param data A data frame containing the input data.
#' @param weights A named vector of weights corresponding to the events.
#' @param conf.int Significance level for calculating confidence intervals. Default is 0.95.
#' @param conf.type Type of confidence interval to calculate. Default is "log".
#' @param testmethod Method to use for the log-rank test. Options are "logrank" or "none". Default is "logrank".
#'
#' @return An object of class "WCEKM" containing the following elements:
#' \describe{
#'   \item{km_data}{A data frame with the combined Kaplan-Meier survival data for all groups.}
#'   \item{events}{A sorted vector of event names where status = 1.}
#'   \item{groups}{A sorted vector of group names.}
#'   \item{sample_size}{A named vector of the sample sizes for each group.}
#'   \item{followup}{The maximum follow-up time across all groups.}
#'   \item{pvalue}{The p-value for the log-rank test.}
#'   \item{testmethod}{The method used for the log-rank test.}
#' }
#' @export
#'
#' @examples
#' # Example usage:
#' weights <- c(CVdeath = 1, MI = 0.55, Stroke = 0.455)
#' WCE_obj <- WCE_KMSurv(surv_exdata, weights)
#'
WCE_KMSurv <- function(data, weights = NULL,
                       conf.int = 0.95,
                       conf.type = "log",
                       testmethod = c( "logrank", "none")) {


  # Data validation
  if (!is.data.frame(data)) {
    stop("The input data must be a data.frame.")
  }

  required_columns <- c("id", "event", "time", "group", "status")
  if (ncol(data) < 5) {
    stop(paste("The input data must have at least 5 columns including:",
               paste(required_columns, collapse = ", ")))
  }

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing in the input data:",
               paste(missing_columns, collapse = ", ")))
  }


  # Check if names of weights match the event names with status = 1
  event_names <- unique(data$event[data$status == 1])
  # Convert 'group' column to ordered
  data$group <- as.ordered(data$group)
  group_lv = levels(data$group)

  # Constant variables
  testmethod <- match.arg(testmethod)
  zval <- qnorm(1 - (1 - conf.int) / 2)
  sample_size <- sapply(split(data$id, data$group), function(x) length(unique(x)))

  # If weights are NULL, set all event weights to 1
  if (is.null(weights)) {
    weights <- setNames(rep(1, length(event_names)), event_names)
    testmethod = "logrank"
  }

  # Process data to calculate weighted events and risk sets
  data <- data %>%
    arrange(time, group, id) %>%
    group_by(group, id) %>%
    mutate(
      event_count = row_number(),
      weight = ifelse(.data$status == 1,
                      weights[as.character(.data$event)] * cumprod(c(1, 1 - head(weights[as.character(.data$event)], -1))), 0),
      n.censor = ifelse(.data$status == 0, 1, 0)
    ) %>%
    ungroup() %>%
    group_by(group) %>%
    mutate(
      n.risk = rev(cumsum(rev(status == 0 | status == 1))),
      n.event = ave(weight, time, FUN = sum)
    ) %>%
    ungroup()%>%
    mutate(n.risk_1 = rev(cumsum(rev(ifelse(group == group_lv[1], 1, 0)))),
           n.risk_2 = rev(cumsum(rev(ifelse(group == group_lv[2], 1, 0)))))


  # Test and calculate p-value
  if (testmethod == "logrank") {
    TestResult = logrankTest(data = data)
    pvalue = TestResult$pvalue
    testmethod = TestResult$method
  } else if (testmethod == "none") {
    TestResult = NULL
    pvalue = NULL
    method = NULL
  } else {
    stop("please input test method!")
  }

  # Generate Kaplan-Meier data
  km_data = survWCEKM(data)

  # If confidence interval type is "log", calculate the confidence intervals for survival probability
  if (conf.type == "log") {
      km_data <- km_data %>%
        mutate(
          lower = exp(log(survival) - zval * std.err), # Calculate lower confidence limit
          upper = pmin(exp(log(survival) + zval * std.err), 1) # Calculate upper confidence limit, limiting to a maximum of 1
        )
  }

  # # Convert group column to original name and order
  km_data$group <- ordered(km_data$group, levels = seq_along(group_lv), labels = group_lv)
  # Create WCEKM object with required elements
  result <- structure(
    list(
      km_data = km_data,
      TestResult = TestResult,
      events = sort(event_names),
      groups = sort(unique(data$group)),
      sample_size = sample_size,
      followup = max(data$time),
      pvalue = pvalue,
      testmethod = testmethod
    ),
    class = "WCEKM"
  )

  return(result)
}
