#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom openxlsx read.xlsx
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @import flexsurv
## usethis namespace: end
NULL

#' Example Survival Data
#'
#' @description
#' `surv_exdata` is an example dataset used for demonstrating survival analysis functions.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{id}{Unique identifier for each subject.}
#'   \item{event}{Event indicator.}
#'   \item{time}{Time to event or censoring.}
#'   \item{group}{Group indicator.}
#'   \item{status}{Status indicator (1 for event, 0 for censored).}
#' }
"surv_exdata"

#' Recurrent Survival Data
#'
#' @description
#' `surv_data_recur` is an example dataset used for demonstrating survival analysis with recurrent events.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{id}{Unique identifier for each subject.}
#'   \item{event}{Event indicator.}
#'   \item{time}{Time to event or censoring.}
#'   \item{group}{Group indicator.}
#'   \item{status}{Status indicator (1 for event, 0 for censored).}
#' }
"surv_data_recur"




#' Example simulation Survival data
#'
#' @description
#' `simulation_data` is an example dataset used for demonstrating survival analysis functions.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{id}{Unique identifier for each subject.}
#'   \item{event}{Event indicator.}
#'   \item{time}{Time to event or censoring.}
#'   \item{group}{Group indicator.}
#'   \item{status}{Status indicator (1 for event, 0 for censored).}
#' }
"simulation_data"
