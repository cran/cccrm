#' Repeated Measures Concordance Correlation Coefficient estimated by U-statistics
#' 
#' @export
#' @description
#' Estimation of the concordance correlation coefficient for repeated measurements using the U-statistics approach. The function is also applicable for the non-repeated measurements scenario.
#' @param dataset An object of class \code{data.frame}.
#' @param ry Character string. Name of the outcome in the data set.
#' @param rmet Character string. Name of the method variable in the data set.
#' @param rtime Character string. Name of the time variable in the data set.
#' @param Dmat  Matrix of weights.
#' @param delta Power of the differences. A value of 0 provides an estimate that is comparable to a repeated measures version of kappa index.
#' @param cl Confidence level.
#' @return A vector that includes the point estimate, confidence interval and standard error of the CCC. Additionally the Fisher's Z-transformation value and its standard error are also provided.
#' @references{
#' King, TS and Chinchilli, VM. (2001). A generalized concordance correlation coefficient for continuous and categorical data. Statistics in Medicine, 20, 2131:2147.
#' 
#' King, TS; Chinchilli, VM; Carrasco, JL. (2007). A repeated measures concordance correlation coefficient. Statistics in Medicine, 26, 3095:3113.
#' 
#' Carrasco, JL; Phillips, BR; Puig-Martinez, J; King, TS;  Chinchilli, VM. (2013). Estimation of the concordance correlation coefficient for repeated measures using SAS and R. Computer Methods and Programs in Biomedicine, 109, 293-304.
#'}
#' @examples
#' # Non-longitudinal scenario
#' newdat=bpres[bpres$NM==1,]
#' estccc=cccUst(newdat,"DIA","METODE")
#' estccc
#' 
# Longitudinal scenario
#' estccc=cccUst(bdaw,"AUC","MET","VNUM")
#' estccc
#' 
# Weighted CCC
#' estccc=cccUst(bfat,"BF","MET","VISITNO",Dmat=diag(c(2,1,1)))
#' estccc
cccUst <-
function(dataset,ry,rmet,rtime=NULL,Dmat=NULL,delta=1,cl=0.95){

# Creating data dataset
dades<-data.frame(dataset)


if (length(rtime)==0) dades <- dades %>% dplyr::rename(y = all_of(ry), 
                                                       met = all_of(rmet))

if (length(rtime)>0) dades <- dades %>% dplyr::rename(y = all_of(ry), 
                                                       met = all_of(rmet),
                                                      time = all_of(rtime))
  
catmet<-unique(dades$met)
if (length(rtime)==0){
  valtime<-"1" 
  ntime=1
  }else {
  valtime<-unique(dades$time)
  ntime<-length(valtime)
  }

ns<-length(dades$y)/(2*ntime)


if (length(Dmat)==0){
  Dmat<-diag(rep(1,ntime))
}

if(sum(dim(Dmat)==c(ntime,ntime))!=2){
stop("Invalid dimensions in weigth matrix")
  }

#Sort data
if (length(rtime)>0){
dades<-dades[order(dades$time),]
}

# Creating Y,X data matrix

Y<-array(dades[dades$met== catmet[1],]$y,c(ns,ntime))
X<-array(dades[dades$met== catmet[2],]$y,c(ns,ntime))

colnames(Y)<-valtime
colnames(X)<-valtime


phimat<-array(NA,c(ns*ns,4))
cont<-0
for (i in 1:ns){
for (j in 1:ns){
cont<-cont+1
phimat[cont,c(3,4)]<-phi(X[i,],Y[i,],X[j,],Y[j,],Dmat,delta) 
phimat[cont,c(1,2)]<-c(i,j)
}
}


colnames(phimat)<-c("i","j","phi1","phi2")

U<-sum((phimat[phimat[,1]!=phimat[,2],3]))/(ns*(ns-1))
V<-sum((phimat[phimat[,1]!=phimat[,2],4]))/(ns*(ns-1))



CCC<-((ns-1)*(V-U))/(U+(ns-1)*V)


phimat1<-phimat[phimat[,1]!=phimat[,2],]
phi1<-tapply(phimat1[,3],phimat1[,1],sum)/(ns-1)
phi2<-tapply(phimat1[,4],phimat1[,1],sum)/(ns-1)

phiv<-cbind(phi1,phi2)
UV<-cbind(U,V)

Saux<-array(0,c(2,2))
C<-array(c(1,0,0,2),c(2,2))
for (i in 1:ns){
Saux<-Saux+t(phiv[i,]-UV)%*%(phiv[i,]-UV)
}

Smat<-C%*%(Saux/(ns^2))%*%C

dev<-array(NA,c(1,2))
dev[1,1]<-((-1)*ns*(ns-1)*V)/((U+(ns-1)*V)^2)
dev[1,2]<-(ns*(ns-1)*U)/((U+(ns-1)*V)^2)

VarCCC<-dev%*%Smat%*%t(dev)

alpha=1-cl

Z<-0.5*(log((1+CCC)/(1-CCC)))
VarZ<-VarCCC/(((1-CCC)^2)*((1+CCC)^2))
ic.z=Z+c(-1,1)*qnorm(1-alpha/2)*c(sqrt(VarZ))
ic.ccc=(exp(2*ic.z)-1)/(exp(2*ic.z)+1)
result<-c(CCC,ic.ccc,sqrt(VarCCC),Z,sqrt(VarZ))
conf.lab=paste((1-alpha)*100,"%",sep="")
names(result)<-c("CCC",paste("LL CI",conf.lab),paste("UL CI",conf.lab),"SE CCC","Z","SE Z")
class(result)<-"cccUst"
result
}

