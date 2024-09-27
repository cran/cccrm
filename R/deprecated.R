

cccvc<-
  function(dataset,ry,rind,rmet,covar=NULL,int=FALSE,cl=0.95,control.lme=list()){
    stopifnot(int %in% c(TRUE,FALSE))
    lifecycle::deprecate_warn("3.0.0", "cccrm::cccvc()", "ccc_vc()",always=TRUE) 
    if (int==TRUE){
          ccc_vc(dataset=dataset,ry=ry,rind=rind,rmet=rmet,int=TRUE,covar=covar,cl=cl,control.lme=control.lme) 
    
    }
    
    else {
      ccc_vc(dataset=dataset,ry=ry,rind=rind,rmet=rmet,int=FALSE,covar=covar,cl=cl,control.lme=control.lme) 
    }
  }




ccclon <-function(dataset,ry,rind,rtime,rmet,covar=NULL,rho=0,cl=0.95,control.lme=list()){
  lifecycle::deprecate_warn("3.0.0", "cccrm::ccclon()", "ccc_vc()",always=TRUE) 
  ccc_vc(dataset=dataset,ry=ry,rind=rind,rtime=rtime,rmet=rmet,covar=covar,rho=rho,cl=cl,control.lme=control.lme) 
}


ccclonw<-function(dataset,ry,rind,rtime,rmet,vecD,covar=NULL,rho=0,cl=0.95,control.lme=list()){
  lifecycle::deprecate_warn("3.0.0", "cccrm::ccclonw()", "ccc_vc()",always=TRUE) 
  ccc_vc(dataset=dataset,ry=ry,rind=rind,rtime=rtime,vecD=vecD,rmet=rmet,covar=covar,rho=rho,cl=cl,control.lme=control.lme) 
}  
  
  
  

