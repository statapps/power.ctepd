\name{power.surv.test}
\alias{power.surv.test}

\title{
Power and sample size for survival distribution (One-sample)
}
\description{
Compute power and sample size of test or determine distribution of alternative hypothesis to obtain target power (same as power.anova.test).
}
\usage{
power.surv.test(n = NULL, sig.level = 0.05, s0 = NULL, s1 = NULL, 
                year = 3, accrual = 3, followup = 3, futile = NULL, 
                futile.prob = 0.01, power = NULL, 
                tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{n}{Number of observations}
  \item{sig.level}{Significance level (Type I error probability)}
  \item{year}{Time point in which survival function is evaluated}
  \item{s0}{Survival function at yaer = year under H0}
  \item{s1}{Survival function at yaer = year under H1}
  \item{accrual}{Years of accrual}
  \item{followup}{Years of followup}
  \item{power}{Power of test (1 minus Type II error probability)}
  \item{futile}{Optinal futility analysis, define the percentage of pearson-year for futility analysis. Value shall be between 0 and 1. Default is NULL.
}
  \item{futile.prob}{Significant level of futility analysis. Probability of rejecting the alternative hypothesis when H1 is true. Default is 0.01}

  \item{tol}{tol for uniroot function}
}

\details{
These calculations use Poisson process and person-year under exponentioal distribution.

Exactly one of the parameters 'n','power' and 's1' must be passed as NULL, and that parameter is determined from the others. 
Notice that the last one has non-NULL default so NULL must be explicitly passed if you want to compute it.
}
\value{
Object of class 'power.htest', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.

  \item{person.year}{Total person-year of followup required for the final analysis}
  \item{median.dur}{Median duration of followup}
  \item{crtl.event}{Critical value for number of observed events to reject the null hypothesis H0}
  \item{crtl.survf}{Critical value for observed survival function at year = year to reject the null hypothesis H0}
  \item{futile.pyear}{Person-year for interim futility analysis}
  \item{futile.event}{Reject H1 if observed event > futile.event}
  \item{futile.prob}{Probability of rejection H1 in the interim analysis when H1 is true}
}
\references{
Lawless, J. Statistical Models and Methods for Lifetime Data, John Wiley and Sons, 1982 (equation 3.2.7 (page 108).
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
\code{\link{power.2surv.test}},
\code{\link{pwr.p.test}}
}
\examples{
##
#pw = power.surv.test(n = 75, s0 = 0.83, s1 = 0.92, sig.level = 0.1)
#pw = power.surv.test(n = 75, s0 = 0.83, power = 0.8, followup = 3.08, sig.level = 0.1)

pw = power.surv.test(power = 0.8, s0 = 0.83, s1 = 0.92, followup = 3, sig.level = 0.1)
print(pw, digits = 4)
}
\keyword{Power} 
\keyword{Sample size}
\keyword{One-sample}
\keyword{Survival}
\keyword{Clinical trial}
\keyword{Phase II}
