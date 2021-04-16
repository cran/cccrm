summary.icc<-
function(object,...){

print(object$model)
cat("\n")
cat("Intraclass Correlation Coefficient - Reliability Index \n")
print(object$icc[1:4])
}


