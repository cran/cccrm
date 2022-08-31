
ic.ccc <-
function(ccc,dev,S,alpha){
  
se.ccc <- as.numeric(sqrt(dev%*%S%*%t(dev)))
z<-ZF(ccc,2)
se.z<-dZF(ccc,2)*se.ccc
ic.z=z+c(-1,1)*qnorm(1-alpha/2)*se.z
ic.icc=(exp(2*ic.z)-1)/(exp(2*ic.z)+1)
result<-c(ccc,ic.icc,se.ccc,z,se.z)
conf.lab=paste((1-alpha)*100,"%",sep="")

names(result)<-c("CCC",paste("LL CI",conf.lab),paste("UL CI",conf.lab),"SE CCC","Z","SE Z")
return(result)
}

