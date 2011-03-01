\name{ccclon}
\alias{ccclon}
\alias{ccclon.default}
\alias{print.ccclon}
\alias{summary.ccclon}
\alias{print.summary.ccclon}
 \title{Longitudinal Repeated Measures Concordance Correlation Coefficient estimated by Variance Components}
\description{Estimates the concordance correlation coefficient for repeated measurements using the variance components from a linear mixed model. }
\usage{
ccclon(dataset,ry,rind,rtime,rmet,covar=NULL,rho=0)
\method{ccclon}{default}(dataset,ry,rind,rtime,rmet,covar=NULL,rho=0)
\method{print}{ccclon}(x,...)
\method{summary}{ccclon}(object,...)
\method{print}{summary.ccclon}(x,...)
}
\arguments{
  \item{dataset}{Name of data set}
  \item{ry}{Character string indicating the outcome in data set }
  \item{rind}{Character string indicating the subject variable in data set}
  \item{rmet}{Character string indicating the method variable in data set}
  \item{rtime}{Character string indicating the time variable in data set} 
  \item{covar}{Character vector indicating the covariables}
  \item{rho}{Within subject correlation structure. The value 0 stands for compound symmetry and 1 for autoregressive order one}  
  \item{x}{Object class ccclon}
  \item{object}{Object class ccclon} 
  \item{\dots}{other arguments to be passed to print or summary}}

\value{
An object of class ccclon.The generic function print and summary gives the estimates of the concordance correlation coefficient.
The object ccclon contains the following components:
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

data(bdaw)

result<-ccclon(bdaw,"cort_auc","SUBJ","VNUM","met")
result

summary(result)
}