\name{cccvc}
\alias{cccvc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Concordance Correlation Coefficient estimated by Variance Components
%%  ~~function to do ... ~~
}
\description{
Estimation of the concordance correlation coefficient for non-repeated measurements and non-longitudinal repeated measurements (replicates) using the variance components from a linear mixed model. The appropriate intraclass correlation coefficient is used as estimator of the concordance correlation coefficient.
}
\usage{
cccvc(dataset, ry, rind, rmet, covar = NULL, int = FALSE, cl=0.95)
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
  \item{rmet}{
Character string. Name of the method variable in the data set.
}
  \item{covar}{
Character vector. Name of covariables to include in the linear mixed model as fixed effects.
}
  \item{int}{
Boolean indicating if the subject-method interaction has to be included in the model
}
  \item{cl}{
Confidence level.
}
}
\details{
The concordance correlation coefficient is estimated using the appropriate intraclass correlation coefficient (see Carrasco and Jover, 2003; Carrasco et al., 2009; Carrasco et al, 2013). The scenarios considered are non-repeated measurements (only one measurement by subject and method) and non-longitudinal repeated measurements, i.e. replicates (multiple measurements by subject and method). The variance components estimates are obtained from a linear mixed model estimated by restricted maximum likelihood. The standard error of CCC is computed using an Taylor's series expansion of 1st order (delta method). Confidence interval is built by applying the Fisher's Z-transformation.
}
\value{
An object of class \code{ccc}. Generic function \code{summary} show a summary of the results. The output is a list with the following components:
\item{ccc}{Concordance Correlation Coefficient estimate}
\item{model}{Summary of the linear mixed model}
\item{vc}{Variance components estimates}
\item{sigma}{Variance components asymptotic covariance matrix}

}

\references{
Carrasco, JL; Jover, L. (2003). Estimating the generalized concordance correlation coefficient through variance components. Biometrics, 59, 849:858.

Carrasco, JL; King, TS; Chinchilli, VM. (2009). The concordance correlation coefficient for repeated measures estimated by variance components. Journal of Biopharmaceutical Statistics, 19, 90:105.

Carrasco, JL; Phillips, BR; Puig-Martinez, J; King, TS;  Chinchilli, VM. (2013). Estimation of the concordance correlation coefficient for repeated measures using SAS and R. Computer Methods and Programs in Biomedicine, 109, 293-304.
}
\author{
Josep Puig-Martinez and Josep L. Carrasco
}



\examples{
data(bpres)

# Scenario 1. 
#Only 1 measure by subject and method. 
#No subject-method interaction included in the model.

newdat=bpres[bpres$NM==1,]
estccc=cccvc(newdat,"DIA","ID","METODE")
estccc
summary(estccc)

# Scenario 2. 
#Two measures by subject and method. 
#No subject-method interaction included in the model.

estccc=cccvc(bpres,"DIA","ID","METODE")
estccc
summary(estccc)

}