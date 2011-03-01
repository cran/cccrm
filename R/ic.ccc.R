ic.ccc <-
function(ccc,dev,S){
se.ccc <- sqrt(dev%*%S%*%t(dev))
z<-0.5*log((1+ccc)/(1-ccc))
se.z<-sqrt(  (se.ccc^2)/(((1+ccc)^2)*((1-ccc)^2))  )
ll.z<-z-qnorm(0.975)*se.z
ul.z<-z+qnorm(0.975)*se.z
ll95<- (exp(2*ll.z)-1)/(exp(2*ll.z)+1)
ul95<- (exp(2*ul.z)-1)/(exp(2*ul.z)+1)
result<-c(ccc,ll95,ul95,se.ccc,z,se.z)
names(result)<-c("CCC","LL CI95%","UL CI95%","SE CCC","Z","SE Z")
return(result)
}

