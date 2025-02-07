resamp_ccc_by_time<-function(i,Bmat,dataset,g_rows,ns,ry,rtime,rmet,vecD, covar, rho, int, cl,control.lme){
  cls <- Bmat[,i]
  xx<-lapply(g_rows[cls],length)
  new_id <- rep(1:ns,times=unlist(xx))
  idx <- unlist(g_rows[cls], recursive = FALSE)
  resamp_data <- dataset[idx, ]
  # Change ind tags
  resamp_data[,"ind"] <- new_id
  resamp_data<-as.data.frame(resamp_data)
  
  ccc_est_by_time(resamp_data, ry, rind="ind", rmet, rtime, plotit = FALSE)$ccc[,2]
  
}