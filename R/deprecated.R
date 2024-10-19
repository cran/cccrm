#' Concordance Correlation Coefficient estimated by Variance Components
#' @export
#' @description
#' Estimation of the concordance correlation coefficient for non-repeated measurements and non-longitudinal repeated measurements (replicates) using the variance components from a linear mixed model. The appropriate intraclass correlation coefficient is used as estimator of the concordance correlation coefficient.
#' @param dataset an object of class \code{data.frame}.
#' @param ry Character string. Name of the outcome in the data set.
#' @param rind Character string. Name of the subject variable in the data set.
#' @param rmet Character string. Name of the method variable in the data set.
#' @param covar Character vector. Name of covariates to include in the linear mixed model as fixed effects.
#' @param int Binary indicating if the subject-method interaction has to be included in the model when analyzing the non-longitudinal setting (defaults to FALSE).
#' @param cl Confidence level.
#' @param control.lme A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help.

#' @details
#' This function has been deprecated. See \link{ccc_vc}.
#' 
cccvc<-
  function(dataset,ry,rind,rmet,covar=NULL,int=FALSE,cl=0.95,control.lme=list()){
    stopifnot(int %in% c(TRUE,FALSE))
    lifecycle::deprecate_warn("3.0.0", "cccrm::cccvc()", "ccc_vc()",always=TRUE) 
    if (int==TRUE){
          ccc_vc(dataset=dataset,ry=ry,rind=rind,rmet=rmet,int=TRUE,covar=covar,cl=cl,control.lme=control.lme) 
    
    }
    
    else {
      ccc_vc(dataset=dataset,ry=ry,rind=rind,rmet=rmet,int=FALSE,covar=covar,cl=cl,control.lme=control.lme) 
    }
  }




#' Concordance Correlation Coefficient estimated by Variance Components
#' @export
#' @description
#' Concordance Correlation Coefficient for longitudinal repeated measures estimated by variance components
#' @param dataset an object of class \code{data.frame}.
#' @param ry Character string. Name of the outcome in the data set.
#' @param rind Character string. Name of the subject variable in the data set.
#' @param rmet Character string. Name of the method variable in the data set.
#' @param rtime Character string. Name of the time variable in the data set.
#' @param covar Character vector. Name of covariates to include in the linear mixed model as fixed effects.
#' @param rho Within subject correlation structure. A value of 0 (default option) stands for compound symmetry and 1 is used for autorregressive of order 1 structure.
#' @param cl Confidence level.
#' @param control.lme A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help.
#' @details
#' This function has been deprecated. See \link{ccc_vc}.
ccclon <-function(dataset,ry,rind,rtime,rmet,covar=NULL,rho=0,cl=0.95,control.lme=list()){
  lifecycle::deprecate_warn("3.0.0", "cccrm::ccclon()", "ccc_vc()",always=TRUE) 
  ccc_vc(dataset=dataset,ry=ry,rind=rind,rtime=rtime,rmet=rmet,covar=covar,rho=rho,cl=cl,control.lme=control.lme) 
}

#' Concordance Correlation Coefficient estimated by Variance Components
#' @export
#' @description
#' Concordance Correlation Coefficient for longitudinal repeated measures estimated by variance components
#' @param dataset an object of class \code{data.frame}.
#' @param ry Character string. Name of the outcome in the data set.
#' @param rind Character string. Name of the subject variable in the data set.
#' @param rmet Character string. Name of the method variable in the data set.
#' @param rtime Character string. Name of the time variable in the data set.
#' @param vecD Vector of weights. The length of the vector must be the same as the number of repeated measures.
#' @param covar Character vector. Name of covariates to include in the linear mixed model as fixed effects.
#' @param rho Within subject correlation structure. A value of 0 (default option) stands for compound symmetry and 1 is used for autorregressive of order 1 structure.
#' @param cl Confidence level.
#' @param control.lme A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help.
#' @details
#' This function has been deprecated. See \link{ccc_vc}.
ccclonw<-function(dataset,ry,rind,rtime,rmet,vecD,covar=NULL,rho=0,cl=0.95,control.lme=list()){
  lifecycle::deprecate_warn("3.0.0", "cccrm::ccclonw()", "ccc_vc()",always=TRUE) 
  ccc_vc(dataset=dataset,ry=ry,rind=rind,rtime=rtime,vecD=vecD,rmet=rmet,covar=covar,rho=rho,cl=cl,control.lme=control.lme) 
}  
  
  
  

