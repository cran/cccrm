icc<-
  function(dataset,ry,rind,covar=NULL,cl=0.95,
           control.lme=list()){
    
    options(contrasts=c("contr.treatment","contr.poly"))
    
    dades<-data.frame(dataset)
    
    dades <- dades %>% dplyr::rename(y = all_of(ry), 
                                     ind = all_of(rind))
    
    dades$ind<-as.factor(dades$ind)
    dades$y<-as.numeric(dades$y)
    
    form=y~1
    if (length(covar)>0){
      form<-as.formula(paste("y~met",paste(covar,sep="+"),sep="+"))}
    
    model.lme<-lme(form,dades,random=~1|ind,method="REML",
                   na.action=na.omit,control=control.lme)
    
    if(is.character(model.lme$apVar)==TRUE){
      stop("Non-positive definite approximate variance-covariance")}
    
    model<-model.lme
    
    # Variance components
    vars<-attr(model$apVar,"Pars")
    SA<-s_exp(vars[1])
    SE<-s_exp(vars[2])
    
     # Calculating icc
    
    icc<-icc0(SA,SE)
    names(icc)<-"ICC"
    
    # Standard error
    alpha=1-cl
    
    D_tau<-matrix(c(d_exp(vars[1]),d_exp(vars[2])),ncol=1)
    
    S<-D_tau%*%t(D_tau)*model$apVar
    D<-matrix(c(d0_1(SA,SE),
                d0_2(SA,SE)),
              nrow=1)
    
    # Number of replicates
    nt<-dades %>% group_by(.data$ind) %>% summarize(n=n())
    ns<-length(unique(nt$ind))
    n<-nrow(dades)
    m<-(n-sum(nt$n^2)/n)/(ns-1)
    
    est<-ic.icc(icc,D,S,0.05,m)
    
    varcomp<-c(SA,SE)
    names(varcomp)<-c("Subject","Random Error")
    
    res<-list(icc=est,vc=varcomp,sigma=S,model=model)
    class(res)<-"icc"
    res 
  }