print.summary.cccUst <-
function(x,...)
{
res<-cbind(x$CCC,x$se,x$low,x$up,x$Z,x$seZ)
colnames(res)<-c("CCC","SE CCC","LL 95%","UL 95%","Z","SE Z")
cat("Results: \n")
print(res)
cat("\n")
}

