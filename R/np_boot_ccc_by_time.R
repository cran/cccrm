#' Non-parametric cluster bootstrap to make inference on the concordance correlation coefficient by time variable
#' @keywords internal
#' @param dataset an object of class \code{data.frame}
#' @param ry Character string. Name of the outcome in the data set.
#' @param rind Character string. Name of the subject variable in the data set.
#' @param rmet Character string. Name of the method variable in the data set.
#' @param rtime Character string. Name of the time variable in the data set.
#' @param covar Character vector. Name of covariates to include in the linear mixed model as fixed effects.
#' @param int Binary indicating if the subject-method interaction has to be included in the model when analyzing the non-longitudinal setting (defaults to FALSE).
#' @param cl Confidence level.
#' @param control.lme A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help.
#' @param nboot Number of bootstrap resamples.
#' @param future_seed Logical/Integer. The seed to be used for parallellization. Further details in \code{furrr_options}.
#' @return Vector of CCC bootstrap estimates.

np_boot_ccc_by_time <- function(dataset,ry, rind, rmet, rtime,
                                covar = NULL, int = FALSE,
                                cl = 0.95, control.lme = list(),
                                nboot = 500,
                                future_seed=123){
  
  
  
  # Resample matrix
  dades<-dataset
  ns <- length(unique(as.character(unlist(dades[,rind]))))
  dots <- as.name(rind)
  grouped <- dplyr::group_by(dades, !!dots)
  g_rows <- dplyr::group_rows(grouped)
  id <- seq_along(g_rows)
  Bmat <- matrix(sample(rep(id,nboot)),nrow=length(id),byrow=FALSE)
  if(anyNA(Bmat)){
    warning("NA values in resample matrix.")
  }
  
  # Compute resamples
  with_progress({
    
    p <- progressr::progressor(steps = nboot)  
    
    ccc_resamples<-furrr::future_map(as.integer(1:nboot), ~ {
      p()
      Sys.sleep(.2)
      resamp_ccc_by_time(.x,Bmat,dades,g_rows,ns,ry,rtime,
                         rmet,vecD, covar, rho, int, cl,control.lme)
    }, .options = furrr::furrr_options(seed = future_seed),
    p = p
    )
    
    
  },enable=TRUE)
  
  k<-length(unique(dades[,rtime]))
  
  ccc_resamples<-matrix(unlist(ccc_resamples),ncol=k,
                        byrow=TRUE)
  
  
  return(ccc_resamples)
}
