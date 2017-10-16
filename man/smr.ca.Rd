\name{smr.ca}
\alias{smr.ca}

\title{
Standard mortality ratio (SMR) for Canadian population
}
\description{
Modify the file from SMR in survexp.fr package to calculation standard mortality ratio for Canadian population.
}
\usage{
smr.ca(time, event, age, sex, ratetable = survexp.ca, alpha = 0.05)
}
\arguments{
  \item{time}{Survival time}
  \item{event}{Event status, 1 for death and 0 for alive}
  \item{age}{Age of the subject}
  \item{sex}{Gender of the subject, 1 for male and 2 for female}
  \item{ratetable}{Expected mortality rate table for Canadian population, default is survexp.ca}
  \item{alpha}{Nominal level of confidence interval. Default is 0.05}
}

\details{
See help file of SMR from 'survexp.fr' for details.
}
\value{
Object of class 'power.htest', a list of the arguments (including the computed one).
  \item{death}{Number of death}
  \item{expected}{Expected number of death in general population}
  \item{SMR.classic}{SMR based on classic method, as in SMR from survexp.fr}
  \item{SMR.poisson}{SMR based on Poisson method, as in SMR from survexp.fr}
  \item{can.haz}{Expected canadian hazard rate of death (matched with age and gender)}
}
\references{
N. Breslow and N. Day, Statistical methods in cancer research, Volume II - The design and analysis of cohort studies, World Health Organization, 1987.
}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}
%%\note{}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{SMR}}, 
}
\examples{
##
# smr = smr.ca(time, event, age, sex)
}
\keyword{smr: standard mortality ratio} 
\keyword{Canada population}
