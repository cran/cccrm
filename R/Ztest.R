#' Wald's test on the Concordance Correlation Coefficient
#'
#' @export
#' @description
#' Estimation of the concordance correlation coefficient for either non-repeated, non-longitudinal, or longitudinal repeated measurements using the variance components from a linear mixed model. The appropriate intraclass correlation coefficient is used as estimator of the concordance correlation coefficient.
#' 
#' @param cccfit An object of class \code{ccc}.
#' @param r0 Integer. Null hypothesis value.
#' @param tr Logical. Should the transformed CCC be used? Only applies if a transformation was applied in the \code{ccc} object. Default value is TRUE.
#' @details
#' A one sided test to the null hypothesis value \deqn{\rho_0}.
#' \deqn{z=\frac{\hat{\theta}-\rho_0}{SE\left(\hat{\theta}\right)}} where \deqn{\hat{\theta}} stands for the CCC estimate and \deqn{SE\left(\hat{\theta}\right)} its standard error.
#' The p-value is computed as \deqn{P\left(X>z\right)} where X follows a standard Normal distribution.
#' @examples
#' 
#' # Testing the CCC is above 0.8
#' ccc_mc=ccc_vc(bpres,"DIA","ID","METODE")
#' ccc_mc
#' Ztest(ccc_mc,r0=0.8)
#' 
#' @return A data frame with two columns: \code{Z}, the statistical test value; and the P-value associated.
Ztest<-function(cccfit,r0=0,tr=TRUE){
  
  stopifnot(tr %in% c(TRUE,FALSE))
  
  if(!is.null(cccfit$resamples)) tr=FALSE
  
  if (tr){
    Z=(cccfit$ccc[5]-atanh(r0))/cccfit$ccc[6]
    pval<-1-pnorm(Z)  
  }else{
    Z=(cccfit$ccc[1]-r0)/cccfit$ccc[4]
    pval<-1-pnorm(Z)
  }
  
  return(data.frame(Z=Z,P.value=pval))
}