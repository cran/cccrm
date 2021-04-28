
Lmat<-function(k){
  aux1<-diag(k)*0
  aux1[1,1:k]<-1
  
  for (i in 1:(k-1)) aux1[1+i,1+i]=1
  L=array(0,dim=c(k,k*(k-1)/2))
  cont=0
  for (i in 1:k){
    for (j in 1:(k-1)){
      if (i>j){
        cont=cont+1
        L[,cont]=aux1[,i]-aux1[,j]
      }
    }
  }	
  return(L)
}

varB<-function(model,nobs) model$varFix[1:nobs,1:nobs]

#Computes variability between observers' means
sb<-function(k,b,varB){
  
  L<-Lmat(k)   # Building L matrix
  
  difmed=t(L)%*%b  # Calculating observers sum of squares
  
  #vardifmed=t(L)%*%varB%*%L
  A<-L%*%t(L)
  aux2<-(t(difmed)%*%difmed)-sum(diag(A%*%varB))
  sb<-max(aux2/(k*(k-1)),0)
  return(sb)
}

# Variance of sb
var.sb<-function(k,b,varB){
  L<-Lmat(k)    # Building L matrix
  A<-L%*%t(L)
  xx<-( 2*sum(diag( ( (A%*%varB)**2 ) ) )+ 4*( t(b)%*%A%*%varB%*%A%*%b ) ) /((k*(k-1))**2)
  return(xx)
}



s_exp<-function(tau)  exp(2*(tau))

d_exp<-Deriv(s_exp,"tau")

icc0<-function(sa,se) sa/(sa+se)
icc1<-function(sa,se,sb) sa/(sa+sb+se)
icc2<-function(sa,sab,se,sb) sa/(sa+sab+sb+se)
icc3<-function(sa,sab,sag,se,sb) (sa+sag)/(sa+sab+sag+se+sb)
icc4<-function(sa,sab,sag,se,sb,sumd) (sumd*(sa+sag))/((sumd*(sa+sab+sag+se))+sb)
  

d0_1<-Deriv(icc0,"sa")
d0_2<-Deriv(icc0,"se")

d1_1<-Deriv(icc1,"sa")
d1_2<-Deriv(icc1,"se")
d1_3<-Deriv(icc1,"sb")

d2_1<-Deriv(icc2,"sa")
d2_2<-Deriv(icc2,"sab")
d2_3<-Deriv(icc2,"se")
d2_4<-Deriv(icc2,"sb")

d3_1<-Deriv(icc3,"sa")
d3_2<-Deriv(icc3,"sab")
d3_3<-Deriv(icc3,"sag")
d3_4<-Deriv(icc3,"se")
d3_5<-Deriv(icc3,"sb")

d4_1<-Deriv(icc4,"sa")
d4_2<-Deriv(icc4,"sab")
d4_3<-Deriv(icc4,"sag")
d4_4<-Deriv(icc4,"se")
d4_5<-Deriv(icc4,"sb")

