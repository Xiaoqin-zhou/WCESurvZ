#' logrankTest
#'
#' Perform a log-rank test for weighted composite endpoint survival analysis
#'
#' This function performs a log-rank test on weighted composite endpoint survival data, calculating the observed and expected events, the chi-square statistic, and the p-value.
#'
#' @param data A data frame containing the input data. The data frame must include the following columns: id, event, time, group, status, weight, n.risk_1, n.risk_2.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{expected}{A vector of expected events for each group.}
#'   \item{observed}{A vector of observed events for each group.}
#'   \item{chi}{The chi-square statistic for the log-rank test.}
#'   \item{pvalue}{The p-value for the log-rank test.}
#' }
#' @export
#'
#' @examples
#' # Example usage:
#' result <- logrankTest(surv_exdata)
#'
logrankTest <- function(data) {

  # Select necessary columns from the input data frame
  data_temp = data %>% select(id, event, time, group, status, weight, n.risk_1, n.risk_2)

  # Initialize variables
  nsize = nrow(data_temp)
  ngroup = length(unique(data_temp$group))
  ord_time = drop(data_temp$time)
  ord_status = drop(data_temp$status)
  ord_group = drop(data_temp$group)
  ord_weight = drop(data_temp$weight)
  ord_n.risk_1 = drop(data_temp$n.risk_1)
  ord_n.risk_2 = drop(data_temp$n.risk_2)
  ord_n.risk_all = ord_n.risk_1 + ord_n.risk_2

  # Call the external C function to perform the log-rank test
  xx <- .C(Cwcelogrank, as.integer(nsize),
           as.integer(ngroup),
           as.double(ord_time),       # time
           as.integer(ord_status),    # status
           as.integer(ord_group),     # group
           as.double(ord_weight),     # weight
           as.double(ord_n.risk_1),   # n.risk_1
           as.double(ord_n.risk_2),   # n.risk_2
           as.double(ord_n.risk_all), # n.risk_all
           observed = double(ngroup),
           expected = double(ngroup),
           tempted = double(ngroup),
           var.e    = double(1),
           double(ngroup))

  # Extract observed and expected events
  otmp = xx$observed
  etmp = xx$expected
  df = (etmp > 0)

  # Calculate chi-square statistic
  temp2 = ((otmp - etmp)[df])[-1]
  vv = xx$var.e
  chi = sum(solve(vv, temp2) * temp2)
  df = (sum(1 * (etmp > 0))) - 1
  pvalue = pchisq(chi, df, lower.tail = FALSE)

  # Return the results as a list
  list(expected = etmp,
       observed = otmp,
       chi = chi,
       pvalue = pvalue)
}
