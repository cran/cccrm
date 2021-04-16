print.icc<-
function(x,...)
{
cat("Intraclass Correlation Coefficient - Reliability Index: \n")
print(x$icc[1:4])
cat("\n")
}
