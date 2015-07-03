\name{cccUst}
\alias{cccUst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Repeated Measures Concordance Correlation Coefficient estimated by U-statistics
%%  ~~function to do ... ~~
}
\description{
Estimation of the concordance correlation coefficient for repeated measurements using the U-statistics approach. The function is also applicable for the non-repeated measurements scenario.
}
\usage{
cccUst(dataset, ry, rmet, rtime = NULL, Dmat = NULL, delta = 1, cl = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
An object of class \code{data.frame}.
}
  \item{ry}{
Character string. Name of the outcome in the data set.
}
  \item{rmet}{
Character string. Name of the method variable in the data set.
}
  \item{rtime}{
Character string. Name of the time variable in the data set.
}
  \item{Dmat}{
Matrix of weights.
}
  \item{delta}{
Power of the differences. A value of 0 provides an estimate that is comparable to a repeated measures version of kappa index.
}
  \item{cl}{
Confidence level.
}
}

\value{
A vector that includes the point estimate, confidence interval and standard error of the CCC. Additionally the Fisher's Z-transformation value and its standard error are also provided.
}
\references{
King, TS and Chinchilli, VM. (2001). A generalized concordance correlation coefficient for continuous and categorical data. Statistics in Medicine, 20, 2131:2147.

King, TS; Chinchilli, VM; Carrasco, JL. (2007). A repeated measures concordance correlation coefficient. Statistics in Medicine, 26, 3095:3113.

Carrasco, JL; Phillips, BR; Puig-Martinez, J; King, TS;  Chinchilli, VM. (2013). Estimation of the concordance correlation coefficient for repeated measures using SAS and R. Computer Methods and Programs in Biomedicine, 109, 293-304.
}
\author{
Josep Puig-Martinez and Josep L. Carrasco
}

\examples{
# Non-longitudinal scenario
data(bpres)
newdat=bpres[bpres$NM==1,]
estccc=cccUst(newdat,"DIA","METODE")
estccc

# Longitudinal scenario
data(bdaw)
estccc=cccUst(bdaw,"AUC","MET","VNUM")
estccc

# Weighted CCC
data(bfat)
estccc=cccUst(bfat,"BF","MET","VISITNO",Dmat=diag(c(2,1,1)))
estccc
}
