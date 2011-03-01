print.summary.cccvc2 <-
function(x,...)
{
cat("Results: \n")
cat("CCC: \n")
print(x$ccc.i)
cat("\n")
cat("Variance Components: \n")
print(x$vc)
cat("\n")
}
