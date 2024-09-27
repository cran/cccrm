resamp_ccc<-function(i,Bmat,dataset,g_rows,ns,ry,rtime,rmet,vecD, covar, rho, int, cl,control.lme){
  cls <- Bmat[,i]
  xx<-lapply(g_rows[cls],length)
  new_id <- rep(1:ns,times=unlist(xx))
  idx <- unlist(g_rows[cls], recursive = FALSE)
  resamp_data <- dataset[idx, ]
  # Change ind tags
  resamp_data[,"ind"] <- new_id
  resamp_data<-as.data.frame(resamp_data)
  
  new_model <- lme_model(resamp_data, ry, rind="ind", rtime, rmet, 
                         vecD, covar, rho, int, cl,control.lme,apVar=FALSE)
  
  #new_model<-try(update(model,data=resamp_data))
  
  
  if(inherits(new_model,"lme")){
    return(ccc_est(new_model, sd_est = FALSE)$ccc)
  }
  
  if(inherits(new_model,"try-error")){
    new_model <- lme_model(resamp_data, ry, rind="ind", rtime, rmet, 
                           vecD, covar, rho, int, cl,
                           control.lme=list(opt = 'optim',minAbsParApVar = 0.1))
    
    
    if(inherits(new_model,"try-error")){
        return(NA)
      }else{
        return(ccc_est(new_model, sd_est = FALSE)$ccc)
      }
    }  
}
