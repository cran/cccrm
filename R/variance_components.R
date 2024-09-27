## Extract relevant variance components of an lme object

variance_components <- function(model){
  vc <- diag(nlme::getVarCov(model))
  vc <- vc |> as.data.frame() |> t() |> as.data.frame() |>
    dplyr::select(!tidyselect::starts_with("met")[-1]) |>
    dplyr::select(!tidyselect::starts_with("time")[-1])
  se <- model$sigma^2
  return(c(unlist(vc),se))
}
