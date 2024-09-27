ic_ccc <-function(ccc,dev,S,alpha, m = 2, transf = "F2", N = NULL){

  if( !(transf %in% c("F","F2","KG","None") ) ){
    message("Options for transformation argument are F,F2,KG or None. 
            Default option F2 is applied.")
  
    transf<-"F2"
    }
  
    se.ccc <- as.numeric(sqrt(dev%*%S%*%t(dev)))

    if(transf == "F2"){
      m <- 2
    }

    if(transf == "F2" | transf == "F"){
      # Fisher Z transform
      z <- Z_F(ccc,m)
      se.z <- dZ_F(ccc,m)*se.ccc
      ic.z <- z+c(-1,1)*qnorm(1-alpha/2)*se.z
      ic.ccc <- (exp(2*ic.z)-1)/(exp(2*ic.z)+m-1)
    }else if(transf == "KG"){
      # Konishi Gupta Z transform
      z <- Z_KG(ccc,m,N)
      bias <-  (7-5*m)/(N*sqrt(18*m*(m-1)))
      se.z <- dZ_KG(ccc,m,N)*se.ccc
      ic.z <- z-bias+c(-1,1)*qnorm(1-alpha/2)*se.z
      ic.ccc <- Z_inv_KG(ic.z,m,N)
    }else{
      ic.ccc <- ccc + c(-1,1)*qnorm(1-alpha/2)*se.ccc
      z <- NA
      se.z <- NA
    }

    result <- c(ccc,ic.ccc,se.ccc,z,se.z)
    conf.lab <- paste((1-alpha)*100,"%",sep="")

    names(result) <- c("CCC",paste("LL CI",conf.lab),
                      paste("UL CI",conf.lab),"SE CCC","Z","SE Z")
    return(result)
}
