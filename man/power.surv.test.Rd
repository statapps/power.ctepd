\name{power.surv.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{power.surv.test}

\title{
Power and sample size calculations for survival distribution (one sample)
}
\description{
Compute power and sample size of test or determine distribution of alternative hypothesis to obtain target power (same as power.anova.test).
}
\usage{
power.surv.test(n = NULL, sig.level = 0.05, s0 = 0.5, s1 = NULL, 
                         year = 3, followup = 3,
                         power = NULL, tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{n}{Number of observations}
  \item{sig.level}{Significance level (Type I error probability)}
  \item{year}{Time point in which survival function is evaluated}
  \item{follow-up}{Years of followup}
  \item{s0}{Survival function at yaer = year under H0}
  \item{s1}{Survival function at yaer = year under H1}
  \item{power}{Power of test (1 minus Type II error probability)}
  \item{tol}{tol for uniroot}
}

\details{
These calculations use Poisson process and person-year under exponentioal distribution.
}
\value{
Object of class '"power.htest"', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.
}
\references{
Lawless, J.  (2014). Survaval Analysis.
}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}
\note{
'uniroot' is used to solve power equation for unknowns, so you may see errors from it, notably about inability to bracket the root when invalid arguments are given.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{power.t.test}}, 
\code{\link{pwr.p.test}}, 
\code{\link{pwr.2p.test}}
}
\examples{
##
#ppt = power.surv.test(n = 75, s0 = 0.83, s1 = 0.92, sig.level = 0.1)
#ppt = power.surv.test(n = 75, s0 = 0.83, power = 0.8, followup = 3.08, sig.level = 0.1)

pn = power.surv.test(power = 0.8, s0 = 0.83, s1 = 0.92, followup = 3, sig.level = 0.1)
print(pn, digits = 4)
}
\keyword{ Power } 
\keyword{ Sample size}
\keyword{ One sample survival distribution}
