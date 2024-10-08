% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccc_sim_data.R
\name{ccc_sim_data}
\alias{ccc_sim_data}
\title{Data simulation using fixed and random effects}
\usage{
ccc_sim_data(
  n = 30,
  nrep = 1,
  nsim = 1,
  model = NULL,
  b = NULL,
  g = NULL,
  mu = 0,
  sa = 1,
  sab = 0,
  sag = 0,
  bg = NULL,
  se = 1,
  future_seed = TRUE,
  workers = 15,
  extra.info = TRUE,
  ...
)
}
\arguments{
\item{n}{Integer. Number of subjects}

\item{nrep}{Integer. Number of replicates}

\item{nsim}{Integer. Number of data sets simulated.}

\item{model}{Object of class \code{lme}.}

\item{b}{Vector. Method fixed effects.}

\item{g}{Vector. Time fixed effects.}

\item{mu}{Integer. Overall mean.}

\item{sa}{Integer. Standard deviation of subject's random effect.}

\item{sab}{Integer. Standard deviation of subject-method interaction's random effect.}

\item{sag}{Integer. Standard deviation of subject-time interaction's random effect.}

\item{bg}{Vector. Method-time interaction's fixed effects. The vector of effects have to be ordered by method and time.}

\item{se}{Integer. Standard deviation of random error effect.}

\item{future_seed}{Logical/Integer. The seed to be used for parallellization. Further details in \code{\link[furrr]{furrr_options}}.}

\item{workers}{Integer. Number of cores to be used for parallellization. Default is 15. Capped to number of available cores minus 1.}

\item{extra.info}{Logical. Should the information about CCC/ICC and variance components simulated be shown?  Default is set to TRUE.}

\item{...}{To pass further arguments.}
}
\value{
A data frame with the simulated data.
}
\description{
The fixed effects and standard deviations of random effects can be set to specific values or, alternatively, obtained from an object of class \code{lme}.
}
\details{
Random effects are simulated as normal distributions with mean 0 and the correspondign standard deviations. The simulated data is obtained
as the addition of the simulated values and the fixed efffects. Parallel computation is used except if data is simulated from an object of class `lme`. In this case.
data is simulated using the \code{\link[nlmeU]{simulateY}} function from \code{nlmeU} package.
}
\examples{
\donttest{
# # Reliability data: 
# 50 subjects, one method, one time, 2 replicates
# Overall mean: -0.25; Subjects standard deviation: 1.5, Random error standard deviation: 1
set.seed(101)
df <- ccc_sim_data(n=50, b = NULL, g = NULL, mu = -0.25, sa = 1.5, se = 1, nrep=2)

# Method comparison data (non-longitudinal)
# 50 subjects, two methods, 2 replicates
# Overall mean: -0.25; Subjects standard deviation: 1.5, Random error standard deviation: 1
# Difference of means between methods 2 and 1: 1
# Three data sets simulated

set.seed(202)
df <- ccc_sim_data(n=50, nsim=3,b = c(0,1), mu = -0.25, sa = 1.5, se = 1, nrep=2)

# Method comparison data (longitudinal)
# 50 subjects, two methods, 3 times, 1 replicate, 
# Overall mean: -0.25; Subjects standard deviation: 1.5, Random error standard deviation: 1
# Difference of means between methods 2 and 1: 1
# Difference of means between times 3,2 and 1 respectively: 0.5 and 0.25.
# Subject-methods interaction standard deviation: 0.25
# Subject-times interaction standard deviation: 0.5
# Same difference of means at each time

set.seed(202)
df <- ccc_sim_data(n=50, b = c(0,1), g=c(0,0.25,0.5), mu = -0.25, sa = 1.5, 
sab=0.25,sag=0.5,se = 1, nrep=2)

# Simulate data using the estimates of a linear mixed model
set.seed(2024)
df3 <- ccc_sim_data(n=50, b = c(0,1), g=c(0,0.25,0.5), mu = -0.25, sa = 1.5, 
                    sab=0.25,sag=0.5,bg=c(0,0.5,0.75,0,1,1),se = 1, nrep=2)
mod3 <- lme_model(df3,"y","id","times","met",control.lme=nlme::lmeControl(opt = 'optim'))
ccc_sim_data(nsim=10,model=mod3)
}

}
\seealso{
\code{\link{ccc_vc}}
}
