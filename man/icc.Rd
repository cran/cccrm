\name{icc}
\alias{icc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Intraclass correlation coefficient
%%  ~~function to do ... ~~
}
\description{
Estimates the one-way intraclass correlation coefficient using the variance components from a linear mixed model. 
}
\usage{
icc(dataset, ry, rind, covar = NULL, cl=0.95, control.lme=list())
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

  \item{covar}{
Character vector. Name of covariables to include in the linear mixed model as fixed effects.
}

  \item{cl}{
Confidence level.
}
  \item{control.lme}{
A list of control values for the estimation algorithm used in \code{lme} function. For further details see \code{lme} help. 
  }
}
\details{
The intraclass correlation coefficient (ICC) (Fleiss, 1986) is estimated using a linear mixed model with subjects (cluster) random effect and random error as components. The effects are assumed to follow Normal distributions with variances Sa and Se respectively. The expression of the ICC is \eqn{Sa/(Sa + Se)}. This intraclass correlation coefficient is also known as reliability index. The variance components estimates are obtained from a linear mixed model (LMM) estimated by restricted maximum likelihood. The function \emph{lme} from package \emph{nlme} (Pinheiro et al., 2021) is used to estimate the LMM. The standard error of CCC is computed using an Taylor's series expansion of 1st order (Ver Hoef, 2012). Confidence interval is built by applying the Fisher's Z-transformation.
}
\value{
An object of class \code{icc}. Generic function \code{summary} show a summary of the results. The output is a list with the following components:
\item{icc}{Intraclass Correlation Coefficient estimate}
\item{model}{nlme object with the fitted linear mixed model}
\item{vc}{Variance components estimates}
\item{sigma}{Variance components asymptotic covariance matrix}

}
\references{



Fleiss, JL. (1986). Reliability of Measurement in *Design and analysis of clinical experiments*. New York:Wiley.

Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2021). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-152, https://CRAN.R-project.org/package=nlme.

Ver Hoef, J.M. (2012) Who Invented the Delta Method?, The American Statistician, 66:2, 124-127.
}

\author{
Josep L. Carrasco
}



\examples{
data(bpres)
 
# Reliability of method 1 

newdat=bpres[(bpres$METODE==1),]
est=icc(newdat,"DIA","ID")
est
summary(est)

# Reliability of method 2 

newdat=bpres[(bpres$METODE==2),]
est=icc(newdat,"DIA","ID")
est
summary(est)


}