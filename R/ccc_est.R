#' Concordance correlation Coefficient estimation from a linear mixed model
#'
#' @export
#' @description It computes the Concordance Correlation Coefficient and its asymptotic confidence interval.
#'
#' @keywords internal
#' @param model The lme model.
#' @param D Weights vector.
#' @param cl Confidence level (0.95 as a default). Bounded between 0 and 1.
#' @param transf Character string. Whether to apply a transformation of the coefficient for inference. Valid options are: "F" for Fisher's Z-transformation; "F2" For Fisher's Z-transformation setting m=2 (default); "KG" Konishi-Gupta transformation; "None", no transformation is applied. See *Details* for further information.
#' @param sd_est Logical. Whether to estimate the asymptotic standard deviation (defaults to TRUE) or to only report the \code{ccc} value.
#' @param ... To pass further arguments.
#'
#' @return A \code{ccc} class object.
#' @seealso \code{\link{ccc_vc}}
#' @examples
#' set.seed(1984)
#' df <- ccc_sim_data(n=50,b = c(0,1), mu = -0.25, sa = 1.5, se = 1, nrep=2)
#' mod <- lme_model(df,"y","id",rmet="met")
#' ccc_est(mod)
ccc_est <- function(model, D = NULL, cl = 0.95, transf = "F2", sd_est = TRUE,...){
  
  # Model data
  dades <- model$data
  
  ns <-length(unique(dades$ind))# Number of subjects
  
  if(!is.null(dades$time)){
    nt <- length(unique(dades$time)) # Number of visits
  }else{
    nt <- 1
  }
  
  if(!is.null(dades$met)){
    nm <- length(unique(dades$met))# Number of observers
    type <- "ccc"
  }else{
    nm <- 1
    type <- "icc"
  }
  
  k <- length(dades$y)/(ns*nt*nm) # Number of replicas
  
  
  
  # Variance components
  vars <- variance_components(model)
  N_end <- length(vars)
  
  SA <-  vars[1]
  
  sumd <- 1
  
  # baseline condition for SA,SAB,SAG,SE,SB
  cond.vc<-c(1,1,1,1,1)
  
  if(type == "icc"){
    cond.vc<-c(1,0,0,0,1)
    SB<-0
  }else if(is.null(dades$time)){
    
    b <- model$coef$fixed[1:nm]
    varB <- model$varFix[1:nm,1:nm]
    SB <- sb(nm,b,varB)
    
    
  }else{
    
    b <- as.matrix(model$coef$fixed)
    
    # Take covariates out
    cond<-c(1:((nm-1)+(nt-1)+1),(length(b)-(nm-1)*(nt-1)+1):length(b))
    b<-b[cond,]
    
    # Design matrix
    
    L <- Lmat_lon(nm,nt,b,dades)
    
    if(!is.null(D)){
      #Weigths matrix
      nd <- nm*(nm-1)/2
      
      if (nd ==1){
        auxD <- D
      }
      
      if (nd > 1) {
        auxD <- bdiag(list(D,D))
        cont <- 2
        while(cont<nd){
          cont <- cont+1
          auxD <- bdiag(list(auxD,D))
        }
      }
      
      sumd <- sum(D)
    }else{
      auxD <- diag(nrow = ncol(L),ncol = ncol(L))
      sumd <- 1
    }
    
    Sb <- model$varFix[cond,cond] # Var-cov of fixed effects
    
    difmed <- t(L)%*%b
    
    A<-L%*%auxD%*%t(L) # A <- L%*%t(L)
    aux1<-(t(difmed)%*%auxD%*%difmed)-sum(diag((A%*%Sb))) #  aux1 <- (t(difmed)%*%difmed)-sum(diag((A%*%Sb)))
    
    SB <- max(aux1/(nm*(nm-1)*nt),0)
  }
  
  if(N_end > 2){
    SAB <- vars[2]
  }else{
    SAB<-0
    cond.vc[2]<-0
  }
  
  if(!is.null(dades$time)){
    SAG <- vars[3]
  }else{
    SAG<-0
    cond.vc[3]<-0
  }
  
  SE <- vars[N_end]
  
  # Variance Components
  
  varcomp<-c(SA,SAB,SAG,SB,SE)
  names(varcomp)<-c("Subjects","Subjects-Method","Subjects-Time","Method","Error")
  varcomp <- varcomp[which(cond.vc==1)]
  
  
  # Calculating icc/ccc
  
  if(sd_est){
    # Standard error
    
    # Model
    
    if( is.character(model$apVar)==TRUE) {
      stop("Approximate variance-covariance matrix not available")
    }
    
    
    if(!is.null(model$modelStruct$corStruct)){
      S_tau <- model$apVar[c(1:3,5),c(1:3,5)] # change for rho = 1
    }else{
      S_tau <- model$apVar
    }
    
    
    if ( ( (cl<0) | (cl>1) ) | (!inherits(cl,"numeric") ) ){
      message("Confidence level must be a numeric value between 0 and 1. 
              Confidence level set to 95%")
      
      cl<-0.95      
    }
    
    alpha <- 1-cl
    if (type=="ccc"){
      
      if(is.null(dades$time)){
        
        var.SB <- var.sb(nm,b,varB)}else{
          
          var.SB <- ((2 * sum(diag(((A %*% Sb)^2)))) + (4 * (t(b) %*% 
                                                               A %*% Sb %*% A %*% b)))/((nm * (nm - 1) * nt)^2)
        }
    }
    
    D_tau <- matrix(2*varcomp[!names(varcomp) %in% "Method"], ncol = 1) # change for rho = 1
    
    # ICC variance estimation
    S <- D_tau%*%t(D_tau)*S_tau
    if(type != "icc"){
      S <- S_mat(S,var.SB,ns,k)
    }
    
    deriv_icc <- c("sa","sab","sag","se","sb")
    derivs <- sapply(deriv_icc, function(i) Deriv::Deriv(icc4,i))
    d_derivs <- sapply(derivs, function(.f) .f(SA,SAB,SAG,SE,SB, sumd))
    
    D <- matrix(d_derivs[which(cond.vc == 1)], nrow = 1)
    
    if(type == "icc"){
      m=k
      icc <- (SA+SAG)/(SA+SAB+SAG+SB+SE)
      names(icc)<-"ICC"
      est <- ic_ccc(icc, D, S, alpha, m, transf, N = ns)
      names(est)[c(1,4)] <- c("ICC","SE ICC")
      res <- list(ccc=est,vc=varcomp,sigma=S,model=model,transf=transf,m=m,N=ns)
      class(res)<-"ccc"
    }else{
      m=nm*nt*k
      ccc <- (SA+SAG)/(SA+SAB+SAG+SB+SE)
      names(ccc)<-"CCC"
      est<-ic_ccc(ccc, D, S, alpha, m, transf, N = ns) #nm*nt*k
      res<-list(ccc=est,vc=varcomp,sigma=S,model=model,transf=transf,m=m,N=ns)
      class(res)<-"ccc"
      
    }
  }else{
    est <- (SA+SAG)/(SA+SAB+SAG+SB+SE)
    res<-list(ccc=est,vc=varcomp)
  }
  
  return(res)
}
