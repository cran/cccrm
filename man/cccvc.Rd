\name{cccvc}
\alias{cccvc}
\alias{cccvc1}
\alias{cccvc1.default}
\alias{print.cccvc1}
\alias{summary.cccvc1}
\alias{print.summary.cccvc1}
\alias{cccvc2}
\alias{cccvc2.default}
\alias{print.cccvc2}
\alias{summary.cccvc2}
\alias{print.summary.cccvc2}
 \title{Concordance Correlation Coefficient by Variance Components
}
\description{Estimates the concordance correlation coefficient through the variance components of a linear mixed model 
}
\usage{
cccvc(dataset,ry,rind,rmet,covar=NULL,int=FALSE)
cccvc1(dataset,ry,rind,rmet,covar=NULL)
\method{cccvc1}{default}(dataset,ry,rind,rmet,covar=NULL)
\method{print}{cccvc1}(x,...)
\method{summary}{cccvc1}(object,...)
\method{print}{summary.cccvc1}(x,...)
cccvc2(dataset,ry,rind,rmet,covar=NULL)
\method{cccvc2}{default}(dataset,ry,rind,rmet,covar=NULL)
\method{print}{cccvc2}(x,...)
\method{summary}{cccvc2}(object,...)
\method{print}{summary.cccvc2}(x,...)
}
\arguments{

  \item{dataset}{Name of data set}
  \item{ry}{Character string indicating the outcome in data set }
  \item{rind}{Character string indicating the subject variable in data set}
  \item{rmet}{Character string indicating the method variable in data set} 
  \item{covar}{Character vector indicating the covariables}
  \item{int}{Boolean value indicating whether the Subject-Method interaction is included in the model} 
  \item{x}{Object class cccvc}
  \item{object}{Object class cccvc} 
  \item{\dots}{other arguments to be passed to print or summary}}
\value{
An object of class cccvc.The generic function print and summary gives the estimates of the concordance correlation coefficient.
The object cccvc contains the following components:
\item{model}{Linear mixed model output}
\item{vc}{Variance Components estimates}
\item{sigma}{An approximate covariance matrix for the variance components}
\item{ccc.p}{The Concordance Correlation Coefficient estimate} 
\item{ccc.i}{Vector containing the Concordance Correlation Coefficient estimate, standard error, 95 percent 
confidence intervals. Additionally Z Fisher's transformation and its standard error are provided}
}
\references{
Carrasco, J. L., Jover, L. (2003). Estimating the generalized concordance correlation coefficient through variance components. Biometrics 59 849 858

Carrasco, J. L., King, T. S., and Chinchilli, V. M. (2009). The concordance correlation coefficient for repeated measures estimated by variance components. Journal of Biopharmaceutical Statistics 19, 90 105.
}
\examples{
data(bpres)

#Case without interaction

result1<-cccvc(bpres,"DIA","ID","METODE")
result1

summary(result1)

#Case adjusted by covariates

result2<-cccvc(bpres,"DIA","ID","METODE",cov=c("EDAD","SEXO"))
result2

summary(result2)

}






