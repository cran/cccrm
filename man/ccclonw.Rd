\name{ccclonw}
\alias{ccclonw}
\alias{ccclonw.default}
\alias{print.ccclonw}
\alias{summary.ccclonw}
\alias{print.summary.ccclonw}
 \title{Weighted Longitudinal Repeated Measures Concordance Correlation Coefficient estimated by Variance Components}
\description{ Estimates the weighted concordance correlation coefficient for repeated measurements using the variance components from a linear mixed model. The different weights are assigned according to the repeated measures.
}
\usage{
ccclonw(dataset,ry,rind,rtime,rmet,vecD,covar=NULL,rho=0)
\method{ccclonw}{default}(dataset,ry,rind,rtime,rmet,vecD,covar=NULL,rho=0)
\method{print}{ccclonw}(x,...)
\method{summary}{ccclonw}(object,...)
\method{print}{summary.ccclonw}(x,...)
 }
\arguments{
  \item{dataset}{Name of data set}
  \item{ry}{Character string indicating the outcome in data set }
  \item{rind}{Character string indicating the subject variable in data set}
  \item{rmet}{Character string indicating the method variable in data set}
  \item{rtime}{Character string indicating the time variable in data set} 
  \item{vecD}{Vector of weigths}
  \item{covar}{Character vector indicating the covariables}
  \item{rho}{Within subject correlation structure. The value 0 stands for compound symmetry and 1 for autoregressive order one}  
  \item{x}{Object class ccclonw}
  \item{object}{Object class ccclonw} 
  \item{\dots}{other arguments to be passed to print or summary}}
\value{
An object of class ccclonw.The generic function print and summary gives the estimates of the concordance correlation coefficient.
The object ccclonw contains the following components:
\item{model}{Linear mixed model output}
\item{vc}{Variance Components estimates}
\item{sigma}{An approximate covariance matrix for the variance components}
\item{ccc.p}{The Concordance Correlation Coefficient estimate} 
\item{ccc.i}{Vector containing the Concordance Correlation Coefficient estimate, standard error, 95 percent 
confidence intervals. Additionally Z Fisher's transformation and its standard error are provided}
}
\references{
King, T. S., Chinchilli, V. M., Carrasco, J. L. (2007). A repeated measures concordance
correlation coefficient. Statistics in Medicine 26(16):3095 3113

Carrasco, J. L., King, T. S., and Chinchilli, V. M. (2009). The concordance correlation coefficient for repeated measures estimated by variance components. Journal of Biopharmaceutical Statistics 19, 90 105.
}
\examples{
data(bfat)

result<-ccclonw(bfat,"bf","SUBJECT","VISITNO","met",vecD=c(2,1,1))
result

summary(result)
}

