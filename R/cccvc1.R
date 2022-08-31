
cccvc1<-
function(dataset,ry,rind,rmet,covar=NULL,cl=0.95,
         control.lme=list()){

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

model.lme<-lme(form,dades,random=~1|ind,method="REML",
               na.action=na.omit,control=control.lme)

if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")}

model<-model.lme

n<-length(unique(dades$ind))# Number of subjects
k<-length(unique(dades$met))# Number of observers
m<-length(dades$y)/(n*k)


# Variance components
vars<-attr(model$apVar,"Pars")
SA<-s_exp(vars[1])
SE<-s_exp(vars[2])

# Calculating observers sum of squares

b<-model$coef$fixed[1:k]
SB<-sb(k,b,varB(model,k))


# Calculating CCC

ccc<-icc1(SA,SE,SB)
names(ccc)<-"CCC"

# Standard error
alpha=1-cl


D_tau<-matrix(c(d_exp(vars[1]),d_exp(vars[2])),ncol=1)

S<-array(NA,c(3,3))
S[1:2,1:2]<-D_tau%*%t(D_tau)*model$apVar
S[3,3]<-var.sb(k,b,varB(model,k))
S[1,3]<-S[3,1]<-(-1/n*m)*S[1,2]
S[2,3]<-S[3,2]<-(-1/n*m)*S[2,2]


D<-matrix(c(d1_1(SA,SE,SB),
            d1_2(SA,SE,SB),
            d1_3(SA,SE,SB)),nrow=1)


est<-ic.ccc(ccc,D,S,0.05)

varcomp<-c(SA,SB,SE)
names(varcomp)<-c("Subject","Observer","Random Error")

res<-list(ccc=est,vc=varcomp,sigma=S,model=model)
class(res)<-"ccc"
res 
}
