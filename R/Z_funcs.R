# Z_KG <- function(r,m,N){
#   sqrt((m-1)/(2*m))*log((1+(m-1)*r)/(1-r)) + (7-5*m)/(N*sqrt(18*m*(m-1)))
# }
# 
# Z_inv_KG <- function(z,m,N){
#   bias <-  (7-5*m)/(N*sqrt(18*m*(m-1)))
#   s <- sqrt((2*m)/(m-1))
#   r <- (exp(s*(z-bias))-1)/(exp(s*(z-bias))+m-1)
#   return(r)
# }

Z_KG <- function(r,m,N){
  sqrt((m-1)/(2*m))*log((1+(m-1)*r)/(1-r)) 
}

Z_inv_KG <- function(z,m,N){
  s <- sqrt((2*m)/(m-1))
  r <- (exp(s*z)-1)/(exp(s*z)+m-1)
  return(r)
}



dZ_KG <- Deriv::Deriv(Z_KG,"r")

Z_F <- function(r,m) 0.5*log((1+(m-1)*r)/(1-r))


dZ_F<- Deriv::Deriv(Z_F,"r")
