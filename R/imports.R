#' @importFrom tidyr pivot_longer
#' @importFrom parallelly availableCores
#' @importFrom progressr progressor with_progress
#' @importFrom purrr map 
#' @importFrom tidyselect starts_with
#' @importFrom nlme lme lmeControl
#' @importFrom stats contrasts formula median na.omit pnorm qnorm quantile rnorm sd time update update.formula
#' @importFrom utils globalVariables
#' @importFrom lifecycle deprecate_warn
#' @importFrom future plan
#' @import dplyr furrr nlmeU Deriv
utils::globalVariables(c("y","ind","met","time"))
NULL
