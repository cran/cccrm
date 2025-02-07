#' @importFrom tidyr pivot_longer
#' @importFrom parallelly availableCores
#' @importFrom progressr progressor with_progress
#' @importFrom purrr map 
#' @importFrom tidyselect starts_with
#' @importFrom nlme lme lmeControl
#' @importFrom stats contrasts formula median na.omit pnorm qnorm quantile rnorm sd time update update.formula var pchisq p.adjust 
#' @importFrom utils globalVariables combn
#' @importFrom lifecycle deprecate_warn
#' @importFrom future plan
#' @import dplyr furrr nlmeU Deriv ggplot2
#' @importFrom MASS ginv
utils::globalVariables(c("y","ind","met","time","CCC","LL95","UL95"))
NULL
