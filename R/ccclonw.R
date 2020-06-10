ccclonw <-
function(dataset,ry,rind,rtime,rmet,vecD,covar=NULL,rho=0,cl=0.95,control.lme=list()){

if (length(vecD) == 0) {
stop("Warning: A vector of weights should be provided")}

dades<-data.frame(dataset)
dades <- dades %>% dplyr::rename(y = all_of(ry), 
                                 ind = all_of(rind), 
                                 met = all_of(rmet),
                                 time = all_of(rtime))
dades$met<-as.factor(dades$met)
dades$time2<-as.numeric(dades$time)
dades$time<-as.factor(dades$time)
dades$y<-as.numeric(dades$y)

form=y~met+time+met*time
if (length(covar)>0){
form<-as.formula(paste("y~met+time+met*time",paste(covar,sep="+"),sep="+"))}


if (length(vecD) != length(unique(dades$time))){
stop("Length of the weight vector must be the number of times")}
D<-diag(vecD)

if ((rho!=0) & (rho!=1)){
stop("Rho must be 0(compound simmetry) or 1 (AR1)")
}

if (rho==0){ 

#Coumpund simmetry model

model.lme<-lme(form,dades,
               random=list(
                 ind=pdBlocked(list(~1,pdIdent(form=~-1+met),pdIdent(form=~-1+time)))),
               control=control.lme)

if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")}
model<-summary(model.lme)

# Variance components

vars<-attr(model$apVar,"Pars")
SA<-s_exp(vars[1])
SAB<-s_exp(vars[2])
SAG<-s_exp(vars[3])
SE<-s_exp(vars[4])

}

if (rho==1){


#AR1 model

model.lme<-lme(form,dades,
               random=list(
                 ind=pdBlocked(list(~1,pdIdent(form=~-1+met),pdIdent(form=~-1+time)))),
               correlation=corAR1(form=~time2|ind/met),control=control.lme)

if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")}
model<-summary(model.lme)

# Variance components
vars<-attr(model$apVar,"Pars")
SA<-s_exp(vars[1])
SAB<-s_exp(vars[2])
SAG<-s_exp(vars[3])
SE<-s_exp(vars[4])
}

# Dimensions
ns<-length(unique(dades$ind))
nt<-length(unique(dades$time))
nm<-length(unique(dades$met))

if ((rho==0) | (rho==1)){

b<-as.matrix(model.lme$coef$fixed)

# Design matrix
nd<-nm*(nm-1)/2
C<-array(0,dim=c(length(b),nt*nm))
k<-0
for (i in 1:nm){
for (j in 1:nt){
k<-k+1
C[,k]<-c(1,contrasts(dades$met)[i,],contrasts(dades$time)[j,],contrasts(dades$met)[i,]*contrasts(dades$time)[j,])
}
}

}


#Weigths matrix


if (nd ==1) auxD<-D

if (nd > 1) {
auxD<-bdiag(list(D,D))
cont<-2
while(cont<nd){
cont<-cont+1
auxD=bdiag(list(auxD,D))
}
}



# Difference between methods matrix
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


alpha=1-cl
Sb<-model.lme$varFix# Var-cov of fixed effects


difmed<-t(L)%*%b
A<-L%*%auxD%*%t(L)

aux1<-(t(difmed)%*%auxD%*%difmed)-sum(diag((A%*%Sb)))
SB<-max(aux1/(nm*(nm-1)),0)
sumd<-sum(D);


# calculating the CCC;
ccc<-icc4(SA,SAB,SAG,SE,SB,sumd)
names(ccc)<-"CCC"

# Variance of between-observers variability;

var.SB<-((2*sum(diag(((A%*%Sb)**2))))+(4*(t(b)%*%A%*%Sb%*%A%*%b)))/((nm*(nm-1))^2)


if (rho == 0){
  STAU <- model$apVar
  D_tau<-matrix(c(d_exp(vars[1]),d_exp(vars[2]),d_exp(vars[3]),d_exp(vars[4])),ncol=1)
}

if (rho == 1){
  STAU <- model$apVar[c(1:3,5),c(1:3,5)]
  D_tau<-matrix(c(d_exp(vars[1]),d_exp(vars[2]),d_exp(vars[3]),d_exp(vars[5])),ncol=1)
}

S<-array(NA,c(5,5))
S[1:4,1:4]<-D_tau%*%t(D_tau)*STAU
S[5,5]<-var.SB
S[1,4]<-S[4,1]<-(-1/ns)*(-1/ns)*(S[1,2]+S[1,3])
S[2,4]<-S[4,2]<-(-1/ns)*(S[2,2]+S[2,3])
S[3,4]<-S[4,3]<-(-1/ns)*(S[3,2]+S[3,3])
S[1,5]<-S[5,1]<-(-1/ns)*(S[1,2]+S[1,4])
S[2,5]<-S[5,2]<-(-1/ns)*(S[2,2]+S[2,4])
S[3,5]<-S[5,3]<-(-1/ns)*(S[3,2]+S[3,4])
S[4,5]<-S[5,4]<-(-1/ns)*(S[4,2]+S[4,4])


D<-matrix(c(d4_1(SA,SAB,SAG,SE,SB,sumd),
            d4_2(SA,SAB,SAG,SE,SB,sumd),
            d4_3(SA,SAB,SAG,SE,SB,sumd),
            d4_4(SA,SAB,SAG,SE,SB,sumd),
            d4_5(SA,SAB,SAG,SE,SB,sumd)),nrow=1)



varcomp<-c(SA,SAB,SAG,SB,SE)
names(varcomp)<-c("Subjects","Subjects-Method","Subjects-Time","Method","Error")
est<-ic.ccc(ccc,D,S,alpha)

res<-list(ccc=est,vc=varcomp,sigma=S,model=model)
class(res)<-"ccc"
res 
}

