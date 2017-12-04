\name{simStudy}
\alias{simStudy}
\title{
Performs a simulation study
}
\description{
Performs a simulation study comparing bias and coverage probability when using either GEE or a semiparametric approach in analyzing accelerometer data
}

\usage{
simStudy(n,numsim,beta,nu,mu,inf,nobs)
}

\arguments{
  \item{n}{Number of individuals in a simulated data.}
  \item{numsim}{Number of simulated datesets.}
  \item{beta}{True coefficient for the binary covariate.}
\item{nu}{Shape and rate parameter for Gamma distribution, in which the subject specific random variable Z_i was generated.}
  \item{mu}{Baseline mean minutes of physical activity per bout.}
    \item{inf}{Whether to generate data with informative observation and censoring times.}
      \item{nobs}{Average number of physical activity bouts.}
}


\value{
A simulated dataset is returned with four columns:
[ID, time, min, x1, phi]. 
}

\author{
Jaejoon Song <jjsong2@mdanderson.org>
}


\keyword{accelerometer}

\examples{
##
## Simulation study when observation and censoring time patterns are noninformative
## Each simulated dataset contains data for 100 individuals
## Two datasets are generated, for illustration purposes
## Expected number of physical activity bouts is 7
##
mysim_ind <- simStudy(n=100,numsim=2,beta=-.4,nu=5,mu=12,inf=FALSE,nobs=7)

##
## Simulation study when observation and censoring time patterns are noninformative
## Each simulated dataset contains data for 100 individuals
## Two datasets are generated, for illustration purposes
## Expected number of physical activity bouts set to 7 when X_i = 1 and Z_i <= 1
## Expected number of physical activity bouts set to 2 when X_i = 0 or Z_i > 1
##
mysim_inf <- simStudy(n=100,numsim=2,beta=-.4,nu=5,mu=12,inf=TRUE,nobs=c(7,2))
}
