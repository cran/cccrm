
cccvc2 <-
function(dataset,ry,rind,rmet,covar=NULL,cl=0.95,control.lme=list()){

options(contrasts=c("contr.treatment","contr.poly"))

dades<-data.frame(dataset)

dades <- dades %>% dplyr::rename(y = all_of(ry), 
                                 ind = all_of(rind), 
                                 met = all_of(rmet))

dades$ind<-as.factor(dades$ind)
dades$met<-as.factor(dades$met)
dades$y<-as.numeric(dades$y)


form=y~met
if (length(covar)>0){
form<-as.formula(paste("y~met",paste(covar,sep="+"),sep="+"))}


model.lme<-lme(form,dades,
               random=list(ind=pdBlocked(list(~1,pdIdent(form=~-1+met)))),
               method="REML",na.action=na.omit,control=control.lme)

if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")}
model<-summary(model.lme)

n<-length(unique(dades$ind))# Number fo subjects
k<-length(unique(dades$met))# Number of observers
m<-length(dades$y)/(n*k)


# Variance components
vars<-attr(model$apVar,"Pars")
SA<-s_exp(vars[1])
SAB<-s_exp(vars[2])
SE<-s_exp(vars[3])

# Calculating observers sum of squares

b<-model$coef$fixed[1:k]
SB<-sb(k,b,varB(model,k))


# Calculating CCC

ccc<-icc2(SA,SAB,SE,SB)
names(ccc)<-"CCC"

# Standard error
alpha=1-cl


D_tau<-matrix(c(d_exp(vars[1]),d_exp(vars[2]),d_exp(vars[3])),ncol=1)

S<-array(NA,c(4,4))
S[1:3,1:3]<-D_tau%*%t(D_tau)*model$apVar
S[4,4]<-var.sb(k,b,varB(model,k))
S[1,4]<-S[4,1]<-(-1/n)*(-1/n)*(S[1,2]+S[1,3])
S[2,4]<-S[4,2]<-(-1/n)*(S[2,2]+S[2,3])
S[3,4]<-S[4,3]<-(-1/n)*(S[3,2]+S[3,3])



D<-matrix(c(d2_1(SA,SAB,SE,SB),
            d2_2(SA,SAB,SE,SB),
            d2_3(SA,SAB,SE,SB),
            d2_4(SA,SAB,SE,SB)),nrow=1)


est<-ic.ccc(ccc,D,S,0.05)

varcomp<-c(SA,SAB,SB,SE)
names(varcomp)<-c("Subject","Subject by observer","Observer","Random Error")


res<-list(ccc=est,vc=varcomp,sigma=S,model=model)
class(res)<-"ccc"
res 

}


