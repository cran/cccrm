#' Power and confidence interval range
#'
#' @export
#' @description Power and confidence interval range obtained by simulation
#' @param n Integer. Number of subjects
#' @param nrep Integer. Number of replicates
#' @param nsim Integer. Number of data sets simulated.
#' @param r0 Integer. Null hypothesis value.
#' @param alpha Type-I error rate.
#' @param model object of class \code{lme}.
#' @param b Vector. Method fixed effects.
#' @param g Vector. Time fixed effects.
#' @param mu Integer. Overall mean.
#' @param sa Integer. Standard deviation of subject's random effect. 
#' @param sab Integer. Standard deviation of subject-method interaction's random effect.
#' @param sag Integer. Standard deviation of subject-time interaction's random effect.
#' @param bg Vector. Method-time interaction's fixed effects
#' @param se Integer. Standard deviation of random error effect.
#' @param extra.info Logical. Should the information about CCC and variance components simulated be shown?  Default is set to TRUE.
#' @param vecD Vector of weights. The length of the vector must be the same as the number of repeated measures.
#' @param covar Character vector. Name of covariates to include in the linear mixed model as fixed effects.
#' @param int Binary indicating if the subject-method interaction has to be included in the model when analyzing the non-longitudinal setting (defaults to FALSE).
#' @param rho Within subject correlation structure. A value of 0 (default option) stands for compound symmetry and 1 is used for autorregressive of order 1 structure.
#' @param cl Confidence level.
#' @param control.lme A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help.
#' @param transf Character string. Whether to apply a transformation of the coefficient for inference. Valid options are: "F" for Fisher's Z-transformation; "F2" For Fisher's Z-transformation setting m=2 (default); "KG" Konishi-Gupta transformasion; "None", no transformation is applied. See *Details* for further information.
#' @param future_seed Logical/Integer. The seed to be used for parallellization. Further details in \code{furrr_options}.
#' @param workers Integer. Number of cores to be used for parallellization. Default is 15. Capped to number of available cores minus 1.
#' @return A data frame with the following components:
#' \itemize{
#' \item \code{n} Number of subjects
#' \item \code{reps} Number of replicates
#' \item \code{CCC}. Median of the CCC estimates.
#' \item \code{Power}. Empirical power computed as proportion of times the null hypothesis is rejected using a type-I error rate of \code{alpha}.
#' \item \code{SEICC}. Average of CCC standard errors.
#' \item \code{SEZ}. Average of transformed CCC standard errors.
#' \item \code{Range IC95}. Average of CCC confidence interval widths.
#' }
#' 
#' @details
#' The power and the range of the confidence interval are computed using the approach suggested in Choudhary and Nagaraja (2018). Data sets are simulated by setting the fixed effects values and the standard deviation of the random effects. The CCC and its standard error are estimated in each data set, along with its 95\% confidence interval and the Wald test \code{\link{Ztest}}.
#' 
#' @references Choudhary, P.K. and Nagaraja, H.N. (2018). Measuring Agreement-Models, Methods, and Applications. John Wiley & Sons
#' @examples
#' \donttest{
#' # Power to test the CCC is above 0.8 with 35 subjects and 4 replicates.
#' # Two methods, three times. Simulated CCC=0.87.
#' sim_pw<-sim_power_ccc(n = 35, nrep=4, nsim=500, r0=0.8, b = c(-0.5,0.5), 
#' g=c(-0.25,0,0.25), mu = -0.25, sa = 4,sab=0.5,sag=1,
#' bg=c(-0.5,-0.25,0.25,-0.5,0.25,0.75),se = 1)
#' }
#' 
#' 
#' 
sim_power_ccc<-function(n = 30, nrep = 2, nsim=300, r0=0,alpha=0.05,model=NULL, b = NULL, g = NULL, mu = 0, sa = 1,
                        sab = 0, sag = 0, bg = NULL, se = 1,  
                        extra.info=TRUE,vecD=NULL,covar=NULL,
                        int=FALSE,rho=0,cl=0.95,control.lme=list(),
                        transf="F2",future_seed=TRUE,workers=15){
  
  if( !(transf %in% c("F","F2","KG","None") ) ){
    message("Options for transformation argument are F,F2,KG or None. 
            Default option F2 is applied.")
    
    transf<-"F2"
  }
  
  if(transf %in% c("F","F2","KG") ){
    tr<-TRUE
  }else{
    tr<-FALSE
  }
  # Parallellization
  ncores <- parallelly::availableCores(omit = 1)
  if(workers >= ncores){
    workers <- ncores
  }
  
    oplan <- future::plan("multisession", workers = workers)
    on.exit(future::plan(oplan))
    message("Parallel execution.")
    message(paste("Number of cores available:",ncores,sep=" "))
    message(paste("Number of cores used:",workers,sep=" "))
    
  
  message("Simulating data...")
  sim_data <- ccc_sim_data(n, nrep, nsim, model=model, b, g, mu, sa, sab, sag, bg, se,
                           future_seed,extra.info=extra.info)
  
  if("met" %in% names(sim_data)) rmet="met" else rmet=NULL
  if("times" %in% names(sim_data)) rtime="times" else rtime=NULL
  
  
  message("Computing estimates...")
  message(paste("Transformation used:",transf))
  progressr::with_progress({
    p <- progressr::progressor(steps = nsim)
    ccc_sim <- furrr::future_map_dfr(1:nsim,~{
      p()
      Sys.sleep(.2)
      df_sam<-sim_data[sim_data$sample==.x,]
      res <- try(ccc_vc(dataset=df_sam,ry="y",rind="id",rtime=rtime,rmet=rmet,vecD,covar,int,rho,cl,
                        control.lme=control.lme, transf,
                        boot=FALSE))
      
      if(inherits(res,"try-error")){
        res <- try(ccc_vc(dataset=df_sam,ry="y",rind="id",rtime=rtime,rmet=rmet,vecD,covar,int,rho,cl,
                          control.lme =
                            nlme::lmeControl(list(opt = 'optim',
                                                  minAbsParApVar = 0.1))))
        if(inherits(res,"try-error")){
          out.test<-as.data.frame(matrix(rep(NA,2),nrow=1))
          out.est<-as.data.frame(matrix(rep(NA,6),nrow=1))
        }
      }else{
        out.test<-Ztest(cccfit=res,r0=r0,info=FALSE)
        out.est<-data.frame(matrix(res$ccc,nrow=1))
      }
      out<-cbind(out.est,out.test)
      return(out)
      
    }, .options = furrr::furrr_options(seed = future_seed),
    p = p)
  })
  
  names(ccc_sim)<-c("CCC","LL","UL","SEICC","Ztrans","SEZ","Ztest","Pval")
  

  exp_CCC<-median(ccc_sim$CCC,na.rm = TRUE)
  exp_power<-mean(ccc_sim$Pval<alpha,na.rm = TRUE)
  exp_SEICC<-sqrt(mean(ccc_sim$SEICC^2,na.rm=TRUE))
  exp_SEZ<-sqrt(mean(ccc_sim$SEZ^2,na.rm=TRUE))
  rangeIC<-ccc_sim$UL-ccc_sim$LL
  exp_Range<-mean(rangeIC,na.rm=TRUE)
  precision=round(qnorm(0.975)*sqrt(exp_power*(1-exp_power)/nsim),2)
  
  paste("Range IC",cl*100,"%",sep="")
  
  result<-data.frame(n,nrep,exp_CCC,exp_power,precision,exp_SEICC,exp_SEZ,exp_Range)
  names(result)<-c("n","reps","CCC","Power","Precision","SEICC","SEZ",
                   paste("Range IC",cl*100,"%",sep=""))
  return(result)
}
