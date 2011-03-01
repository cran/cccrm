\name{cccUst}
\alias{cccUst}
\alias{cccUst.default}
\alias{print.cccUst}
\alias{summary.cccUst}
\alias{print.summary.cccUst}
 \title{Repeated Measures Concordance Correlation Coefficient estimated by U statistics}
\description{ Estimates the concordance correlation coefficient for repeated measurements using the U statistics approach.  This function is also appropiate for the non-repeated measures setting subject and observer.
}
\usage{
cccUst(dataset,ry,rmet,rtime,Dmat=NULL,delta=1)
\method{cccUst}{default}(dataset,ry,rmet,rtime,Dmat=NULL,delta=1)
\method{print}{cccUst}(x,...)
\method{summary}{cccUst}(object,...)
\method{print}{summary.cccUst}(x,...)
 }
\arguments{
  \item{dataset}{Name of data set}
  \item{ry}{Character string indicating the outcome in data set }
  \item{rmet}{Character string indicating the method variable in data set}
  \item{rtime}{Character string indicating the time variable in data set} 
  \item{Dmat}{Matrix of weigths}
  \item{delta}{delta}
  \item{x}{CCC results}
  \item{object}{CCC results}
  \item{\dots}{other arguments to be passed to print or summary}}

\value{
An object of class cccUst.The generic function print and summary gives the estimates of the concordance correlation coefficient.
The object cccUst contains the following components:
\item{CCC}{Coficient of concordance}
\item{se}{Standard Error of CCC}
\item{low}{Lower limit of the confidence interval ccc}
\item{up}{upper limit of the confidence interval ccc} 
\item{Z}{Z trasformation}
\item{seZ}{Standard Error of Z trasformation}}
\references{

King, T. S. and Chinchilli, V. M. (2001). A generalized concordance correlation  coefficient for continuous and categorical data. Statistics in Medicine 20, 2131 2147.

King, T. S., Chinchilli, V. M., Carrasco, J. L. (2007). A repeated measures concordance
correlation coefficient. Statistics in Medicine 26(16):3095 3113
}
\examples{
#No Longitudinal dataset
data(bpres)
bp1<-bpres[bpres$NM==1,]

result<-cccUst(bp1,"DIA","METODE")
result

summary(result)

#Longitudinal dataset
data(bdaw)

result2<-cccUst(bdaw,"cort_auc","met","VNUM")
result2

summary(result2)

#Weigthed lognitudinal dataset
data(bfat)

result3<-cccUst(bfat,"bf","met","VISITNO",Dmat=diag(c(2,1,1)))
result3

summary(result3)
}


