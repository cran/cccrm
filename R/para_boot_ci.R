#'  Parametric bootstrap to make inference on the concordance correlation coefficient
#' @keywords internal
#' @param dataset an object of class \code{data.frame}
#' @param model an object of class \code{lme}
#' @param ry Character string. Name of the outcome in the data set.
#' @param rind Character string. Name of the subject variable in the data set.
#' @param rmet Character string. Name of the method variable in the data set.
#' @param rtime Character string. Name of the time variable in the data set.
#' @param vecD Vector of weights. The length of the vector must be the same as the number of repeated measures.
#' @param covar Character vector. Name of covariates to include in the linear mixed model as fixed effects.
#' @param int Binary indicating if the subject-method interaction has to be included in the model when analyzing the non-longitudinal setting (defaults to FALSE).
#' @param rho Within subject correlation structure. A value of 0 (default option) stands for compound symmetry and 1 is used for autorregressive of order 1 structure.
#' @param cl Confidence level.
#' @param control.lme A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help.
#' @param nboot Number of bootstrap resamples.
#' @param parallel Whether the code is parallellized. The parallellization method is \code{multisession}.
#' @param future_seed Logical/Integer. The seed to be used for parallellization. Further details in \code{furrr_options}.
#' @return Vector of CCC bootstrap estimates.
para_boot_ci <- function(dataset,model, ry, rind,
                         rtime = NULL, rmet = NULL, vecD = NULL,
                         covar = NULL, rho = 0, int = F, cl = 0.95,
                         control.lme = list(),
                         nboot = 500, parallel=TRUE,future_seed = 121){
  
  # Generate R simulations of the response values
  ry_sim <- nlmeU::simulateY(model, nsim = nboot, verbose = F)
  row.names(ry_sim) <- 1:model$dims$N
  
  # Pre-calc
  dat <- dataset
  if(!is.null(vecD)){
    D <- diag(vecD)
  }else{
    D <- NULL
  }
  
  # Generate the bootstrap replicates
  with_progress({
    p <- progressr::progressor(steps = nboot)
    
    if(parallel){
      
      ccc_resamples <- furrr::future_map(1:nboot,~{
        p()
        Sys.sleep(.2)
        dat[,ry] <- ry_sim[,.x]
        MOD <- try(lme_model(dat,ry,rind, rtime, rmet, vecD, covar, rho,int, cl,
                             control.lme,apVar=FALSE))
        
        if(inherits(MOD,"try-error")){
          MOD <- try(lme_model(dat,ry,rind, rtime, rmet, vecD, covar, rho,int, cl,
                               control.lme =
                                 nlme::lmeControl(opt = 'optim',
                                                  minAbsParApVar = 0.1),apVar=FALSE))
          if(inherits(MOD,"try-error")){
            return(NA)
          }
        }
        
        res <- try(ccc_est(MOD, D, cl, sd_est = FALSE)$ccc)
        if(inherits(res,"try-error")){
          return(NA)
        }
        return(res)
      }, .options = furrr::furrr_options(seed = future_seed),
      p = p)
      
    }
    
    if(!parallel){
      
      ccc_resamples <- purrr::map(1:nboot,~{
        p()
        Sys.sleep(.2)
        dat[,ry] <- ry_sim[,.x]
        MOD <- try(lme_model(dat,ry,rind, rtime, rmet, vecD, covar, rho,int, cl,
                             control.lme,apVar=FALSE))
        
        if(inherits(MOD,"try-error")){
          MOD <- try(lme_model(dat,ry,rind, rtime, rmet, vecD, covar, rho,int, cl,
                               control.lme =
                                 nlme::lmeControl(opt = 'optim',
                                                  minAbsParApVar = 0.1),apVar=FALSE))
          if(inherits(MOD,"try-error")){
            return(NA)
          }
        }
        
        res <- try(ccc_est(MOD, D, cl, sd_est = FALSE)$ccc)
        if(inherits(res,"try-error")){
          return(NA)
        }
        return(res)
      }, .progress=TRUE)
    }
    
  })
  
  return(unlist(ccc_resamples))
}
