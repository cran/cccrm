cccvc<-
function(dataset,ry,rind,rmet,covar=NULL,int=FALSE,cl=0.95,control.lme=list()){
if (int==TRUE)cccvc2(dataset,ry,rind,rmet,covar,cl,control.lme=list()) 
  else cccvc1(dataset,ry,rind,rmet,covar,cl,control.lme=list())}