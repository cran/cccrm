\name{ccclon}
\alias{ccclon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Concordance Correlation Coefficient for longitudinal repeated measures estimated by variance components}
%%  ~~function to do ... ~~

\description{Estimation of the concordance correlation coefficient for repeated measurements using the variance components from a linear mixed model. The appropriate intraclass correlation coefficient is used as estimator of the concordance correlation coefficient.}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

\usage{
ccclon(dataset, ry, rind, rtime, rmet, covar = NULL, rho = 0, cl = 0.95,
control.lme=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
An object of class \code{data.frame}.
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
The concordance correlation coefficient is estimated using the appropriate intraclass correlation coefficient (see Carrasco et al, 2009; Carrasco et al, 2013). The variance components estimates are obtained from a linear mixed model estimated by restricted maximum likelihood. The standard error of CCC is computed using an Taylor's series expansion of 1st order (delta method). Confidence interval is built by applying the Fisher's Z-transformation.
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
}
\author{
Josep Puig-Martinez and Josep L. Carrasco
}



\seealso{
\code{\link{ccclonw}}
}
\examples{
data(bdaw)
estccc<-ccclon(bdaw,"AUC","SUBJ","VNUM","MET")
estccc
summary(estccc)

}

