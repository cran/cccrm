ccclon <-
function(dataset,ry,rind,rtime,rmet,covar=NULL,rho=0,cl=0.95,
         control.lme=list()){

dades<-data.frame(dataset)

dades <- dades %>% dplyr::rename(y = all_of(ry), 
                                 ind = all_of(rind), 
                                 met = all_of(rmet),
                                 time = all_of(rtime))

dades$ind<-as.factor(dades$ind)
dades$met<-as.factor(dades$met)
dades$y<-as.numeric(dades$y)
dades$time2<-as.numeric(dades$time)
dades$time<-as.factor(dades$time)

form=y~met+time+met*time
if (length(covar)>0){
form<-as.formula(paste("y~met+time+met*time",paste(covar,sep="+"),sep="+"))}


if ((rho!=0) & (rho!=1)){
  stop("Rho must be 0(compound simmetry) or 1 (AR1)")
}


if (rho==0){ 

#Compound simmetry model

model.lme<-lme(form,dades,
               random=list(ind=pdBlocked(list(~1,pdIdent(form=~-1+met),
                                              pdIdent(form=~-1+time)))),
               control=control.lme,na.action=na.omit)

if(is.character(model.lme$apVar)==TRUE){
  stop("Non-positive definite approximate variance-covariance")}
model<-model.lme

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
                 ind=pdBlocked(
                   list(~1,pdIdent(form=~-1+met),pdIdent(form=~-1+time)))),
               correlation=corAR1(form=~time2|ind/met),
               control=control.lme,na.action=na.omit)


if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")
  }
model<-model.lme

# Variance components

# Variance components

vars<-attr(model$apVar,"Pars")
SA<-s_exp(vars[1])
SAB<-s_exp(vars[2])
SAG<-s_exp(vars[3])
SE<-s_exp(vars[5])
}

# Dimensions
ns<-length(unique(dades$ind))
nt<-length(unique(dades$time))
nm<-length(unique(dades$met))

if ((rho==0) | (rho==1)){

b<-as.matrix(model.lme$coef$fixed)

# Take covariates out
cond<-c(1:((nm-1)+(nt-1)+1),(length(b)-(nm-1)*(nt-1)+1):length(b))  
b<-b[cond,]

# Design matrix

L<-Lmat_lon(nm,nt,b,dades)

Sb<-model.lme$varFix[cond,cond]# Var-cov of fixed effects
difmed<-t(L)%*%b
A<-L%*%t(L)

aux1<-(t(difmed)%*%difmed)-sum(diag((A%*%Sb)))
SB<-max(aux1/(nm*(nm-1)*nt),0)
}


# calculating the CCC;

ccc<-icc3(SA,SAB,SAG,SE,SB)
names(ccc)<-c("CCC")


# Variance of between-observers variability;

var.SB<-((2*sum(diag(((A%*%Sb)**2))))+(4*(t(b)%*%A%*%Sb%*%A%*%b)))/((nm*(nm-1)*nt)^2)


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
S[1,5]<-S[5,1]<-(-1/ns)*(S[1,2]+S[1,4])
S[2,5]<-S[5,2]<-(-1/ns)*(S[2,2]+S[2,4])
S[3,5]<-S[5,3]<-(-1/ns)*(S[3,2]+S[3,4])
S[4,5]<-S[5,4]<-(-1/ns)*(S[4,2]+S[4,4])



# dev: Vector of derivatives;
alpha=1-cl

D<-matrix(c(d3_1(SA,SAB,SAG,SE,SB),
            d3_2(SA,SAB,SAG,SE,SB),
            d3_3(SA,SAB,SAG,SE,SB),
            d3_4(SA,SAB,SAG,SE,SB),
            d3_5(SA,SAB,SAG,SE,SB)),nrow=1)



varcomp<-c(SA,SAB,SAG,SB,SE)
names(varcomp)<-c("Subjects","Subjects-Method","Subjects-Time","Method","Error")
est<-ic.ccc(ccc,D,S,alpha)

res<-list(ccc=est,vc=varcomp,sigma=S,model=model)
class(res)<-"ccc"
res 

}