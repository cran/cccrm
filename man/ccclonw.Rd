\name{ccclonw}
\alias{ccclonw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Concordance Correlation Coefficient for longitudinal repeated measures estimated by variance components
%%  ~~function to do ... ~~
}
\description{
Estimation of the concordance correlation coefficient (CCC) for repeated measurements using the variance components from a linear mixed model. The appropriate intraclass correlation coefficient is used as estimator of the concordance correlation coefficient. Weights are assigned to repeated measurements in the CCC computation process.
}
\usage{
ccclonw(dataset, ry, rind, rtime, rmet, vecD, covar = NULL, rho = 0, cl = 0.95, 
control.lme=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
an object of class \code{data.frame}.
}
  \item{ry}{
Character string. Name of the outcome in the data set.
}
  \item{rind}{
Character string. Name of the subject variable in the data set.
}
  \item{rtime}{
Character string. Name of the time variable in the data set.
}
  \item{rmet}{
Character string. Name of the method variable in the data set.
}
  \item{vecD}{
Vector of weigths. The length of the vector must be the same as the number of repeated measures.
}
    \item{covar}{
Character vector. Name of covariables to include in the linear mixed model as fixed effects.
}
  \item{rho}{
Within subject correlation structure. A value of 0 (default option) stands for compound simmetry and 1 is used for autoregressive of order 1 structure.
}
  \item{cl}{
Confidence level.
}
  \item{control.lme}{
A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help. 
  }
}


\details{
The concordance correlation coefficient is estimated using the appropriate intraclass correlation coefficient  which expression is modified accordingly to assign different weights to each repeated measurement (see Carrasco et al, 2009; Carrasco et al, 2013). The variance components estimates are obtained from a linear mixed model (LMM) estimated by restricted maximum likelihood. The function \emph{lme} from package \emph{nlme} (Pinheiro et al., 2021) is used to estimate the LMM. The standard error of CCC is computed using an Taylor's series expansion of 1st order (Ver Hoef, 2012). Confidence interval is built by applying the Fisher's Z-transformation.
}
\value{
An object of class \code{ccc}. Generic function \code{summary} show a summary of the results. The output is a list with the following components:
\item{ccc}{Concordance Correlation Coefficient estimate}
\item{model}{Summary of the linear mixed model}
\item{vc}{Variance components estimates}
\item{sigma}{Variance components asymptotic covariance matrix}
}
\references{
Carrasco, JL; King, TS; Chinchilli, VM. (2009). The concordance correlation coefficient for repeated measures estimated by variance components. Journal of Biopharmaceutical Statistics, 19, 90:105.

Carrasco, JL; Phillips, BR; Puig-Martinez, J; King, TS;  Chinchilli, VM. (2013). Estimation of the concordance correlation coefficient for repeated measures using SAS and R. Computer Methods and Programs in Biomedicine, 109, 293-304.

Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2021). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-152, https://CRAN.R-project.org/package=nlme.

Ver Hoef, J.M. (2012) Who Invented the Delta Method?, The American Statistician, 66:2, 124-127.
}
\author{
Josep Puig-Martinez and Josep L. Carrasco
}


\seealso{
\code{\link{ccclon}}
}
\examples{
data(bfat)
estccc<-ccclonw(bfat,"BF","SUBJECT","VISITNO","MET",vecD=c(2,1,1))
estccc
summary(estccc)
}

