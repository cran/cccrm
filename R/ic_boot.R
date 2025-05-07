ic_boot <- function(resamples,cl = 0.95,boot_ci = "BCa",orig_ccc = NULL){
  # Compute statistics
  alpha <- 1 - cl
  ccc_resamples <- unlist(resamples)
  
  n_NA<-sum(is.na(ccc_resamples))
  
  if (n_NA>0){
    message(paste(n_NA,"resamples removed because NAs"))
    ccc_resamples<-ccc_resamples[!is.na(ccc_resamples)]    
  }

  

  if((boot_ci == "BCa")){

    ci_resamples <- bca(ccc_resamples,orig_ccc,cl)

  }else{

    delta_resamples <- ccc_resamples - orig_ccc
    ci_resamples <- orig_ccc - quantile(delta_resamples,
                                        c(1-alpha/2,alpha/2))
  }

  response <- c(orig_ccc,ci_resamples,
                sd(ccc_resamples, na.rm = T),mean(ccc_resamples),
                median(ccc_resamples, na.rm = T))

  conf.lab <- paste(cl*100,"%",sep="")
  names(response) <- c("CCC",paste("LL CI",conf.lab),
                       paste("UL CI",conf.lab),"SE CCC","mean CCC","median CCC")
  return(response)
}
