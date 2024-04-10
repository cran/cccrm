
ic.icc <-
function(icc,dev,S,alpha,m){
  
se.ccc <- as.numeric(sqrt(dev%*%S%*%t(dev)))
z<-ZF(icc,m)
se.z<-dZF(icc,m)*se.ccc
ic.z=z+c(-1,1)*qnorm(1-alpha/2)*se.z

ic.icc<-inv.ZF(ic.z,m)

result<-c(icc,ic.icc,se.ccc,z,se.z)
conf.lab=paste((1-alpha)*100,"%",sep="")

names(result)<-c("ICC",paste("LL CI",conf.lab),paste("UL CI",conf.lab),"SE ICC","Z","SE Z")
return(result)
}

