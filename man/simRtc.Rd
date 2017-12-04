\name{simRtc}
\alias{simRtc}
\title{
Simulates minutes of physical activity in continuos time
}
\description{
Simulates minutes of physical activity in continuous time
}

\usage{
simRtc(n,beta,nu,mu,inf,nobs)
}

\arguments{
  \item{n}{Number of individuals in the simulated data.}
  \item{beta}{True coefficient for the binary covariate.}
  \item{nu}{Shape and rate parameter for Gamma distribution, in which the subject specific random variable Z_i was generated.}
  \item{mu}{Baseline mean minutes of physical activity per bout.}
  \item{inf}{Whether to generate data with informative observation and censoring times.}
      \item{nobs}{Average number of physical activity bouts.}
}


\value{
A simulated dataset is returned with four columns:
[ID, time, min, x1, phi]. 
The number of physical activity minutes 'min' is generated as mu*Z_i*exp(X_i*beta), where the expected number of bouts is 'nobs'. 
}

\author{
Jaejoon Song <jjsong2@mdanderson.org>
}


\keyword{accelerometer}

\examples{
##
## Simulauting data for a 50 individuals 
## with noninformative observation and censoring time patterns
##
simdata <- simRtc(n=50,beta=-.4,nu=5,mu=12,inf=FALSE,nobs=7)

##
## Simulauting data for a 50 individuals 
## with informative observation and censoring time patterns
##
simdata <- simRtc(n=50,beta=-.4,nu=5,mu=12,inf=TRUE,nobs=c(7,2))
}
