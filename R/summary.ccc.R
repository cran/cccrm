#' @export
summary.ccc<-
function(object,...){

print(object$vc)
cat("\n")
cat("CCC estimated by variance compoments \n")
print(object$ccc[1:4])
}


