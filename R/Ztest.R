#' Wald's test on the Concordance Correlation Coefficient
#'
#' @export
#' @description
#' Wald's test is applied to assess whether the CCC (ICC) is greater than a reference value.
#' Additionally, Wald's test is also used to compare two independent CCC (ICC).
#' 
#' @param cccfit An object of class \code{ccc}.
#' @param cccfit2 An object of class \code{ccc}.
#' @param r0 Integer. Null hypothesis value.
#' @param info Logical. Should information about the transformation used be printed?
#' @details
#' If only one ccc is provided, the function runs a one sided test to the null hypothesis value \deqn{\rho_0}.
#' \deqn{z=\frac{\hat{\theta}-\rho_0}{SE\left(\hat{\theta}\right)}} where \deqn{\hat{\theta}} stands for the CCC estimate and \deqn{SE\left(\hat{\theta}\right)} its standard error.
#' If a second CCC is provided, the function runs a two-sided test to the null hypothesis of equality of CCCs.
#' \deqn{z=\frac{\hat{\theta_1}-\hat{\theta_2}}{\sqrt{Var\left(\hat{\theta_1}\right)}+Var\left(\hat{\theta_2}\right)}}.
#' In both cases, the p-value is computed as \deqn{P\left(X>z\right)} where X follows a standard Normal distribution.
#' The test uses the transformation indicated when the \code{ccc} object was generated. 
#'
#' 
#' @examples
#' 
#' # Testing the CCC is above 0.8
#' ccc_mc=ccc_vc(bpres,"DIA","ID","METODE")
#' ccc_mc
#' Ztest(ccc_mc,r0=0.8)
#' 
#' # Comparing two CCC
#' 
#' bpres_Male <- bpres |> dplyr::filter(SEXO==1)
#' bpres_Female <- bpres |> dplyr::filter(SEXO==2)
#' 
#' ccc_DIA_Male=ccc_vc(bpres_Male,"DIA","ID","METODE")
#' ccc_DIA_Female=ccc_vc(bpres_Female,"DIA","ID","METODE")
#' Ztest(ccc_DIA_Male,ccc_DIA_Female)
#' 
#' @return A data frame with two columns: \code{Z}, the statistical test value; and the P-value associated.
Ztest<-function(cccfit,cccfit2=NULL,r0=0,info=TRUE){
  
if(is.null(cccfit2)){
  transf<-cccfit$transf
  stopifnot(transf %in% c("F","F2","KG","None"))
  m<-cccfit$m
  if(transf == "F2") m<-2
  N<-cccfit$N
  
  if(transf == "F2" | transf == "F"){
    Z=(cccfit$ccc[5]-Z_F(r0,m))/cccfit$ccc[6]
    pval<-1-pnorm(Z)
  }else if(transf == "KG"){
    bias <-  (7-5*m)/(N*sqrt(18*m*(m-1)))
    Z_rho<-Z_KG(r0,m,N)-bias
    Z=(cccfit$ccc[5]-Z_rho)/cccfit$ccc[6]
    pval<-1-pnorm(Z)
  }  
  else{
    Z=(cccfit$ccc[1]-r0)/cccfit$ccc[4]
    pval<-1-pnorm(Z)
  }
  if (info) message(paste("Transformation used:",transf))
  return(data.frame(Z=Z,P.value=pval))
}else{
  
  transf<-cccfit$transf
  transf2<-cccfit2$transf
  stopifnot(transf %in% c("F","F2","KG","None"))
  stopifnot(transf2 %in% c("F","F2","KG","None"))
  mr=(cccfit$ccc[1]+cccfit2$ccc[1])/2
  if (transf!=transf2){
    message("Different transformations used. Test is run with no transformation applied")
    transf<-transf2<-"None"
  }
  m1<-cccfit$m
  m2<-cccfit2$m
  N1<-cccfit$N
  N2<-cccfit2$N
  
  if(transf == "F2") m1<-m2<-2
  
  if(transf == "F2" | transf == "F"){
    bias<-0.5*log((1+(m1-1)*mr)/(1+(m2-1)*mr))
    Z=(cccfit$ccc[5]-cccfit2$ccc[5]-bias)/sqrt(cccfit$ccc[6]^2+cccfit2$ccc[6]^2)
    pval<-1-pnorm(abs(Z)  )
  }else if(transf == "KG"){
    
    bias1 <-  (7-5*m1)/(N1*sqrt(18*m1*(m1-1)))
    bias2 <-  (7-5*m2)/(N2*sqrt(18*m2*(m2-1)))
    
    bias <- sqrt((m1-1)/(2*m1))*log((1+(m1-1)*mr)/(1-mr))-
      sqrt((m2-1)/(2*m2))*log((1+(m2-1)*mr)/(1-mr))+bias1-bias2
    
    
    Z=(cccfit$ccc[5]-cccfit2$ccc[5]-bias)/sqrt(cccfit$ccc[6]^2+cccfit2$ccc[6]^2)
    pval<-1-pnorm(abs(Z))
  }  
  else{
    Z=(cccfit$ccc[1]-cccfit2$ccc[1])/sqrt(cccfit$ccc[4]^2+cccfit2$ccc[4]^2)
    pval<-1-pnorm(abs(Z))
  }
  message(paste("Transformation used:",transf))
  return(data.frame(Z=abs(Z),P.value=pval))
  
  
}

}
