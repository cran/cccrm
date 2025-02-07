#' Concordance Correlation Coefficient estimation by variance components.
#'
#' @export
#' @description
#' Estimation of the non-longitudinal concordance correlation coefficient at each time using the variance components approach.
#' @param dataset an object of class \code{data.frame}.
#' @param ry Character string. Name of the outcome in the data set.
#' @param rind Character string. Name of the subject variable in the data set.
#' @param rmet Character string. Name of the method variable in the data set.
#' @param rtime Character string. Name of the time variable in the data set.
#' @param covar Character vector. Name of covariates to include in the linear mixed model as fixed effects.
#' @param int Binary indicating if the subject-method interaction has to be included in the model when analyzing the non-longitudinal setting (defaults to FALSE).
#' @param cl Confidence level.
#' @param control.lme A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help.
#' @param transf Character string. Whether to apply a transformation of the coefficient for inference. Valid options are: "F" for Fisher's Z-transformation; "F2" For Fisher's Z-transformation setting m=2 (default); "KG" Konishi-Gupta transformation; "None", no transformation is applied. See *Details* for further information.
#' @param future_seed Logical/Integer. The seed to be used for parallellization. Further details in  \code{\link[furrr]{furrr_options}}.
#' @param workers Integer. Number of cores to be used for parallellization. Default is 15. Capped to number of available cores minus 1.
#' @param plotit Logical. If TRUE it generates a plot with the CCC and their confidence intervals for each time. 
#' @param test Logical. If TRUE the equality of CCCs is assessed. Default to FALSE.
#' @param nboot Number of bootstrap resamples.
#' @param adj.method Character string. Correction method for pairwise comparisons. See \code{\link[stats]{p.adjust}}
#' @param ... To pass further arguments.
#'
#' @return A \code{ccc} class object. Generic function \code{summary} show a summary of the results. The output is a list with the following components:
#' \itemize{
#'  \item \code{ccc}. CCC estimates at each level of time variable.
#'  \item \code{plot}. Plot of the CCC along with their confidence intervals.
#'  \item \code{res_test}. Test of equality of CCCs.
#'  \item \code{ph_table}. Pairwise comparison of CCCs.
#'  }
#'  
#' @references{
#'  
#' Davison A.C., Hinkley D.V. (1997). Bootstrap Methods and Their Application. Cambridge: Cambridge University Press.
#'  
#' Field, C.A., Welsh, A.H. (2007). Bootstrapping Clustered Data. Journal of the Royal Statistical Society. Series B (Statistical Methodology) 69(3):369-390.
#' 
#' Vanbelle S. (2017). Comparing dependent kappa coefficients obtained on multilevel data. Biometrical Journal 59(5):1016-1034. 
#'  
#' }
#' @details
#' The concordance correlation coefficient is estimated using the variance components approach. Confidence intervals are built using the asymptotic Normal distribution approach.
#' Variance-covariance matrix of CCC estimates is estimated by non-parametric balanced randomized cluster bootstrap approach (Davison and Hinkley, 1997; Field and Welsh, 2007).
#' Overall equality of CCCs is tested following the non-parametric bootstrap approach suggested in  Vanbelle (2017).
#' 
#' 
#' 
#' @examples
#' 
#' \dontrun{
#' ccc_est_by_time(bdaw, "AUC", "SUBJ", "MET", "VNUM")
#' ccc_est_by_time(bpres, "SIS", "ID", "METODE", "NM",test=TRUE)
#' }
ccc_est_by_time<-function(dataset,ry,rind,rmet, rtime,
                          covar=NULL,int = F,cl=0.95,
                          control.lme=list(), future_seed=TRUE,
                          transf = "F2",workers = 15,plotit=TRUE,test=FALSE,
                          nboot=500,adj.method="holm",
                          ...){
  
  
  dataset[,"time"]<-factor(dataset[,rtime])
  data_sp<-split(dataset,dataset$time)  
  
  out<-sapply(1:length(data_sp),function(i){
    ccc_vc(data_sp[[i]], ry=ry, rind=rind, rmet=rmet,int=int,transf=transf)$ccc[1:3]
  }
  )
  
  results<-data.frame(time=unique(dataset$time),t(out))
  names(results)<-c(rtime,"CCC","LL95","UL95")
  
  if(plotit){
    
    lim1<-max(0,min(results$LL95)-0.1)
    lim2<-min(max(results$UL95)+0.1,1)
    
    plot <- ggplot(results, aes(x=.data[[rtime]], y = CCC)) +
      geom_point() +
      geom_errorbar(aes(ymin = LL95, ymax = UL95), width = 0.2) +
      labs(title = "CCC with 95% Confidence Intervals",
           x = quo_name(rtime),
           y = "CCC") + ylim(lim1,lim2)
    
    theme_minimal() 
    
  } else plot=NULL
  
  if(test){
    # Parallellization
    ncores <- parallelly::availableCores(omit = 1)
    if(workers >= ncores){
      workers <- ncores
      oplan <- future::plan("multisession", workers = workers)
      on.exit(future::plan(oplan))
      message("Parallel execution.")
      message(paste("Number of cores available:",ncores,sep=" "))
      message(paste("Number of cores used:",workers,sep=" "))
      
    }
    
    ccc_resamples<-np_boot_ccc_by_time(dataset,ry,rind,rmet,rtime,covar,int,
                            cl, control.lme,nboot,
                            future_seed)
    
    b<-results$CCC
    S<-var(ccc_resamples)
    k <- length(b)
    
    # Matrix of pairwise differences
    M <- matrix(0, nrow = choose(k, 2), ncol = k)
    row_index <- 1
    
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
        M[row_index, i] <- 1
        M[row_index, j] <- -1
        row_index <- row_index + 1
      }
    }
    
    # Pairwise differences
    L<-M %*% b
    P<-M%*%S%*%t(M)
    cc<-t(L)%*%MASS::ginv(P)%*%L
    df<-k-1
    pval<-1-pchisq(cc,df)
    
    out_test<-data.frame(Chi.Sq=cc,DF=df,Pvalue=pval)
    
    pvals<-1-pchisq((L^2)/diag(P),1)
    pvals_adj<-p.adjust(pvals,method=adj.method)
    
    
    V<-as.character(results[,1])
    combinations <- combn(V, 2, FUN = function(x) paste(x, collapse = "-"))
    labs <- as.vector(combinations)
    
    ph_table<-data.frame(Difs=labs,Estimate=L,SE=sqrt(diag(P)),Adj.P=pvals_adj)
    
    
    
    
  } else {
    out_test=NULL; ph_table=NULL}
  
  return(list(ccc = results, plot = plot, 
              res_test=out_test,pair_comp=ph_table))
}
