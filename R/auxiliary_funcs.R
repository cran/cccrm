
## Building contrast matrices for generic and longitudinal cases

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

  # L matrix for ccclon
  Lmat_lon<-function(nm,nt,b,dades){
    # Design matrix

    # Number of between-methods differences
    nd<-nm*(nm-1)/2
    # All differences design matrix
    C<-array(0,dim=c(length(b),nt*nm))
    k<-0
    for (i in 1:nm){
      for (j in 1:nt){
        k<-k+1
        C[,k]<-c(1,contrasts(dades$met)[i,],contrasts(dades$time)[j,],
                 c(contrasts(dades$met)[i,]%*%t(contrasts(dades$time)[j,])))
      }
    }

    # Difference between methods at each time matrix
    L<-array(0,dim=c(length(b),nt*nd))
    k<-0
    for (i in 1:(nt*(nm-1))){
      for (j in (1:(nm-1))){
        if ((i+nt*j)<=(nt*nm)){
          k<-k+1
          L[,k]=C[,i]-C[,i+nt*j]
        }
      }
    }
    return(L)
  }

# For derivatives

 icc4<-function(sa,sab,sag,se,sb,sumd) (sumd*(sa+sag))/((sumd*(sa+sab+sag+se))+sb)


# Smat_computation

 S_mat <- function(S,var.SB,ns,k){
   init_row <- 2
   if(nrow(S) == 3){
     init_row <- 3
   }
   S_add <- sapply(S[,2], function(x) (-1/ns*k)*x)
   if(nrow(S) == 4){
     S_add <- sapply(S[,2]+S[,4], function(x) (-1/ns*1)*x)
   }
   S_amp <- S |> rbind(S_add)
   S_add <- c(S_add,var.SB)
   S_amp <- S_amp |> cbind(S_add)

   return(S_amp)
 }

# SB auxiliaries

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

 
 phi <-
   function(X1,Y1,X2,Y2,Dmat,delta){
     
     if (delta!=0){
       
       phi1<-0.5*(((abs(X1-Y1))**delta)%*%Dmat%*%((abs(X1-Y1))**delta)+((abs(X2-Y2))**delta)%*%Dmat%*%((abs(X2-Y2))**delta))
       phi2<-0.5*(((abs(X1-Y2))**delta)%*%Dmat%*%((abs(X1-Y2))**delta)+((abs(X2-Y1))**delta)%*%Dmat%*%((abs(X2-Y1))**delta))
       
       return(c(phi1,phi2))
     }
     if (delta==0){
       phi1<-0.5*(((abs(X1-Y1))!=0)%*%Dmat%*%((abs(X1-Y1))!=0)+((abs(X2-Y2))!=0)%*%Dmat%*%((abs(X2-Y2))!=0))
       phi2<-0.5*(((abs(X1-Y2))!=0)%*%Dmat%*%((abs(X1-Y2))!=0)+((abs(X2-Y1))!=0)%*%Dmat%*%((abs(X2-Y1))!=0))
       
       return(c(phi1,phi2))
       
       
     }
   }
 
 bdiag <-
   function(x){
     if(!is.list(x)) stop("x not a list")
     n <- length(x)
     if(n==0) return(NULL)
     x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
       stop("Zero-length component in x"))
     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,]
     cc <- d[2,]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1,-1] <- rcum[-n]
     ind[2,] <- rcum
     ind[3,-1] <- ccum[-n]
     ind[4,] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                                  (y[3]+1):y[4]], imat=imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     return(out)
   }
 
 
 bca<-function (theta, t0, cl = 0.95)
 {
   
   theta<-theta[!is.na(theta)]
   cl_low <- (1 - cl)/2
   cl_hi <- 1 - cl_low
   nsims <- length(theta)
   z.inv <- length(theta[theta < t0])/nsims
   z <- qnorm(z.inv)
   U <- (nsims - 1) * (t0 - theta)
   A1 <- sum(U^3)
   A2 <- 6 * (sum(U^2))^{
     3/2
   }
   a <- A1/A2
   ll.inv <- pnorm(z + (z + qnorm(cl_low))/(1 - a * (z + qnorm(cl_low))))
   ll <- quantile(theta, ll.inv, names = FALSE)
   ul.inv <- pnorm(z + (z + qnorm(cl_hi))/(1 - a * (z + qnorm(cl_hi))))
   ul <- quantile(theta, ul.inv, names = FALSE)
   return(c(ll, ul))
 }
