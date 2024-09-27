#' Fits a Linear Mixed Effects Model
#'
#' @param dataset an object of class \code{data.frame}.
#' @param ry Character string. Name of the outcome in the data set.
#' @param rind Character string. Name of the subject variable in the data set.
#' @param rtime Character string. Name of the time variable in the data set.
#' @param rmet Character string. Name of the method variable in the data set.
#' @param vecD Vector of weights. The length of the vector must be the same as the number of repeated measures.
#' @param covar Character vector. Name of covariates to include in the linear mixed model as fixed effects.
#' @param rho Within subject correlation structure. A value of 0 (default option) stands for compound symmetry and 1 is used for autoregressive of order 1 structure.
#' @param int Boolean indicating if the subject-method interaction has to be included in the model.
#' @param cl Confidence level.
#' @param control.lme A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help.
#' @param apVar Logical. Should the asymptotic variance-covariance matrix of the variance components be estimated in the linear mixed model? (Defaults to TRUE). 
#' @param ... To pass further arguments.
#'
#' @return an object of class \code{lme}.
#' @export
#'
#' @examples
#' # Reliability ICC
#' set.seed(2024)
#' df <- ccc_sim_data(b = NULL, g = NULL, mu = -0.25, sa = 1.5, se = 1)
#' mod1 <- lme_model(df,"y","id")
#' mod1
#'
#' #Non-longitudinal Methods comparison data
#' set.seed(2024)
#' df2 <- ccc_sim_data(n=50,b = c(0,1), mu = -0.25, sa = 1.5, se = 1, nrep=2)
#' mod2 <- lme_model(df2,"y","id",rmet="met")
#' mod2
#' 
#' # Longitudinal Methods comparison data
#' set.seed(2024)
#' df3 <- ccc_sim_data(n=50, b = c(0,1), g=c(0,0.25,0.5), mu = -0.25, sa = 1.5, 
#'                     sab=0.25,sag=0.5,bg=c(0,0.5,0.75,0,1,1),se = 1, nrep=2)
#' 
#' mod3 <- lme_model(df3,"y","id","times","met",control.lme=nlme::lmeControl(opt = 'optim'))
#' mod3
#'
lme_model <- function(dataset,ry,rind,rtime = NULL,rmet = NULL,vecD = NULL,
                      covar=NULL,rho=0,int = FALSE,cl=0.95,control.lme=list(),
                      apVar=TRUE,...){
  
  stopifnot(rho %in% c(0,1))
  stopifnot(is.logical(int))
  
  dades<- data.frame(dataset) |> dplyr::select(any_of(c(ry,rind,rtime,rmet,covar)))
  
  dades <- dades |> dplyr::rename(y = all_of(ry),
                                  ind = all_of(rind),
                                  met = all_of(rmet),
                                  time = all_of(rtime))
  
  rand_form <- ~1|ind
  cor_form <- NULL
  
  dades$y<- as.numeric(dades$y)
  dades$ind<- as.factor(dades$ind)
  
  if(!is.null(rmet)){
    dades$met <- as.factor(dades$met)
    form <- y~met
    
    if(int){
      rand_form <- list(ind=nlme::pdBlocked(list(~1,
                                                 nlme::pdIdent(form=~-1+met))))
    }
    
  }else{
    form <- y~1
  }
  
  if(!is.null(rtime)){
    dades$time2 <- as.numeric(dades$time)
    dades$time <- as.factor(dades$time)
    
    if(!is.null(rmet)){
      form <- update.formula(form,~.+time+met*time)
    }else{
      form <- update.formula(form,~.+time)
    }
    
    rand_form <- nlme::reStruct(
      list(ind=nlme::pdBlocked(list(~1,nlme::pdIdent(form=~-1+met),
                                    nlme::pdIdent(form=~-1+time)))),
      data = dades)
    
    if(rho == 1){
      dades <- dades |> group_by(ind,met,time) |> mutate(k=seq_along(time))
      
      if (sum(dades$k>1)>0){
        message("For AR1 option time covariate must have unique values within groups. 
                Argument rho is set to 0")
        
      }else{
        cor_form <- nlme::corAR1(form=~time2|ind/met)
      }
    }
  }else{
    options(contrasts=c("contr.treatment","contr.poly"))
  }
  
  if(!is.null(covar)){
    covar_update <- paste("~ . +",paste(covar,collapse =" + "))
    
    form <- update(form,covar_update)
  }
  
  if (!apVar) control.lme$apVar<-FALSE
  
  model.lme <- nlme::lme(form,dades,random= rand_form, correlation = cor_form,
                         method="REML",na.action=na.omit,control=control.lme)
  
  model.lme$call$fixed <- form
  model.lme$call$random <- formula(model.lme$modelStruct$reStruct)
  model.lme$call$correlation <- cor_form
  model <- model.lme
  
  if( (apVar) & (is.character(model.lme$apVar)==TRUE) ){
    stop("Non-positive definite approximate variance-covariance")
  }
  
  return(model)
}
