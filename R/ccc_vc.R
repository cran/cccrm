#' Concordance Correlation Coefficient estimation by variance components.
#'
#' @export
#' @description
#' Estimation of the concordance correlation coefficient for either non-repeated, non-longitudinal, or longitudinal repeated measurements using the variance components from a linear mixed model. The appropriate intraclass correlation coefficient is used as estimator of the concordance correlation coefficient.
#' 
#' @param dataset an object of class \code{data.frame}.
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
#' @param transf Character string. Whether to apply a transformation of the coefficient for inference. Valid options are: "F" for Fisher's Z-transformation; "F2" For Fisher's Z-transformation setting m=2 (default); "KG" Konishi-Gupta transformation; "None", no transformation is applied. See *Details* for further information.
#' @param boot Logical. Whether to compute the CCC confidence interval by bootstrapping or asymptotic methods (defaults to FALSE).
#' @param boot_param Logical. Whether to compute a parametric bootstrap or a non-parametric bootstrap (defaults to FALSE).
#' @param boot_ci Character. Type of bootstrap confidence interval. Either "BCa" (which is the default) or "empirical".
#' @param nboot Integer. Number of bootstrap resamples. Default is 300.
#' @param parallel Logical. Whether the code is parallellized. The parallellization method is \code{multisession}.
#' @param future_seed Logical/Integer. The seed to be used for parallellization. Further details in  \code{\link[furrr]{furrr_options}}.
#' @param workers Integer. Number of cores to be used for parallellization. Default is 15. Capped to number of available cores minus 1.
#' @param sd_est Logical. Whether to estimate the asymptotic standard deviation (defaults to TRUE) or to only report the \code{ccc/icc} value.
#' @param apVar Logical. Should the asymptotic variance-covariance matrix of the variance components be estimated in the linear mixed model? (Defaults to TRUE). 
#' @param ... To pass further arguments.
#'
#' @return A \code{ccc} class object. Generic function \code{summary} show a summary of the results. The output is a list with the following components:
#' \itemize{
#'  \item \code{ccc}. CCC/ICC estimate
#'  \item \code{model}. nlme object with the fitted linear mixed model.
#'  \item \code{vc}. Variance components estimates.
#'  \item \code{sigma}. Variance components asymptotic covariance matrix.
#'  }
#' @details
#' The concordance correlation coefficient is estimated using the appropriate intraclass correlation coefficient (see Carrasco and Jover, 2003; Carrasco et al., 2009; Carrasco et al, 2013). 
#' 
#' The scenarios considered are: a) reliability assessment (several measurements taken with one method); b) methods comparison data with non-repeated measurements (only one measurement by subject and method);
#' c) Methods comparison data with non-longitudinal repeated measurements, i.e. replicates (multiple measurements by subject and method); and d) Methods comparison data with longitudinal repeated measurements (multiple longitudinal measurements by subject and method).
#' 
#'  The variance components estimates are obtained from a linear mixed model (LMM) estimated by restricted maximum likelihood. The function \emph{lme} from package \emph{nlme} (Pinheiro et al., 2021) is used to estimate the LMM. 
#'  
#'  The standard error of CCC and its confidence interval can be obtained: a) asymptotically, using Taylor's series expansion of 1st order (Ver Hoef, 2012); b) using balanced randomized cluster bootstrap approach (Davison and Hinkley, 1997; Field and Welsh, 2007); c) using parametric bootstrap (Davison and Hinkley, 1997).
#'  
#'  When estimating asymptotically the standard error, the confidence intervals are built using the point estimate of the CCC/ICC, its standard error, and the appropriate quantile of the  standard Normal distribution.
#'  However, the approximation to the asymptotic Normal distribution is improved if the CCC/ICC is transformed using the Fisher's Z-transformation (Fisher, 1925), or the Konishi-Gupta transformation (Konishi and Gupta, 1989).
#'  In case the number of replicates is equal to 2, both transformations give the same result.
#'  
#' 
#' @references{
#' 
#' Carrasco, JL; Jover, L. (2003). Estimating the generalized concordance correlation coefficient through variance components. Biometrics, 59, 849:858.
#'
#' Carrasco, JL; King, TS; Chinchilli, VM. (2009). The concordance correlation coefficient for repeated measures estimated by variance components. Journal of Biopharmaceutical Statistics, 19, 90:105.
#'
#' Davison A.C., Hinkley D.V. (1997). Bootstrap Methods and Their Application. Cambridge: Cambridge University Press.
#'
#' Field, C.A., Welsh, A.H. (2007). Bootstrapping Clustered Data. Journal of the Royal Statistical Society. Series B (Statistical Methodology). 69(3), 369-390.
#'
#' Fisher, R. A. (1925) Statistical Methods for Research Workers. Edinburgh: Oliver
#'
#' Konishi, S. and Gupta, A. K. (1989) Testing the equality of several intraclass correlation coefficients. J Statist. Planng Inf., 21, 93-105.
#'
#'
#' Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2021). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-152, \url{https://CRAN.R-project.org/package=nlme}.
#'
#' Ver Hoef, J.M. (2012) Who Invented the Delta Method?, The American Statistician, 66:2, 124-127.
#' } 
#' @examples
#' 
#' \dontrun{
#' # Scenario 1. Reliability 
#' newdat <- bpres |> dplyr::filter(METODE==1)
#' icc_rel<-ccc_vc(newdat,"DIA","ID")
#' icc_rel
#' summary(icc_rel)
#' 
#' 
#' # Confidence interval using non-parametric bootstrap
#'
#' icc_rel_bt<-ccc_vc(newdat,"DIA","ID",boot=TRUE,sd_est=FALSE,
#' nboot=500,parallel=TRUE)
#' icc_rel_bt
#' summary(icc_rel_bt)
#' 
#' 
#' #' # Scenario 2. Non-longitudinal methods comparison.
#' # Only 1 measure by subject and method. 
#' # No subjects-method interaction included in the model.
#' 
#' newdat <- bpres |> dplyr::filter(NM==1)
#' ccc_mc<-ccc_vc(newdat,"DIA","ID","METODE")
#' ccc_mc
#' summary(ccc_mc)
#' 
#' 
#' # Confidence interval using parametric bootstrap
#'
#' ccc_mc_bt<-ccc_vc(newdat,"DIA","ID",boot=TRUE,boot_param=TRUE,
#' sd_est=FALSE,nboot=500,parallel=TRUE)
#' ccc_mc_bt
#' summary(ccc_mc_bt)
#' 
#
#' 
#' # Scenario 3. Non-longitudinal methods comparison. 
#' # Two measures by subject and method. 
#' # No subject-method interaction included in the model.
#' 
#' ccc_mc=ccc_vc(bpres,"DIA","ID","METODE")
#' ccc_mc
#' summary(ccc_mc)
#' 
#' 
#' 
#' # Scenario 4. Methods comparison in longitudinal repeated measures setting.
#' ccc_mc_lon<-ccc_vc(bdaw,"AUC","SUBJ","MET","VNUM")
#' ccc_mc_lon
#' summary(ccc_mc_lon)
#' 
#' # Scenario 5. Methods comparison in longitudinal repeated measures setting.
#' # More weight given to readings from first time.
#' 
#' ccc_mc_lonw<-ccc_vc(bfat,"BF","SUBJECT","MET","VISITNO",vecD=c(2,1,1))
#' ccc_mc_lonw
#' summary(ccc_mc_lonw)
#' }
#' 
ccc_vc <- function(dataset,ry,rind,rmet = NULL,rtime = NULL,vecD = NULL,
                      covar=NULL,int = F,rho=0,cl=0.95,
                      control.lme=list(), transf = "F2",
                      boot = FALSE, boot_param = FALSE,
                      boot_ci = "BCa", nboot = 300, parallel = FALSE,future_seed = TRUE,
                      workers = 15, sd_est=TRUE,apVar=TRUE,...){

  if(!is.null(vecD)){
    if (length(vecD) != length(unique(dataset[[rtime]]))){
      stop("Length of the weight vector must be the number of times")
    }

    D<-diag(vecD)

  }else{
    D <- NULL
  }

  model <- lme_model(dataset, ry, rind, rtime, rmet, vecD, covar, rho, int, cl,
                     control.lme,apVar)

  # Parallellization
  ncores <- parallelly::availableCores(omit = 1)
  if(workers >= ncores){
    workers <- ncores
  }

  if(parallel){
    oplan <- future::plan("multisession", workers = workers)
    on.exit(future::plan(oplan))
    message("Parallel execution.")
    message(paste("Number of cores available:",ncores,sep=" "))
    message(paste("Number of cores used:",workers,sep=" "))
  }

  # CCC statistic + Confidence Interval
  if(boot){

    # Bootstrap resamples
    if(boot_param){
      resamples <- para_boot_ci(dataset,model, ry, rind, rtime, rmet, vecD,
                               covar, rho, int, cl,control.lme,nboot, future_seed,parallel = parallel)
    }else{
      resamples <- np_boot_ci(dataset, ry, rind, rtime, rmet, vecD, covar, rho, int, cl,
                             control.lme,nboot)
    }
    names(resamples)<-NULL
    # Bootstrap CI
    orig_ccc <- try(ccc_est(model, D, cl,transf, sd_est = F,control.lme=control.lme))
    if(inherits(orig_ccc,"try-error")){
      stop("Failed to compute CCC estimate.")
    }
      
    res <- ic_boot(resamples,cl,boot_ci,orig_ccc$ccc)
    response <-list(ccc=res,vc=orig_ccc$vc,sigma=NULL,model=model,resamples=resamples,
                    transf=NULL,m=NULL,N=NULL)
    class(response)<-"ccc"

  }else{
    # Asymptotic
    response <- ccc_est(model, D, cl,transf=transf,sd_est=sd_est)
  }


  return(response)
}

