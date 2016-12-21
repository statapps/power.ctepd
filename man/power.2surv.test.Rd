\name{power.2surv.test}
\alias{power.2surv.test}

\title{
Power and sample size for survival distribution (Two-sample)
}
\description{
Compute power and sample size of test or determine distribution of alternative hypothesis to obtain target power (same as power.surv.test).
}
\usage{
power.2surv.test(n = NULL, sig.level = 0.05, hr = NULL, s0 = 0.7,
          year = 3, accrual = 3, followup = 3, power = NULL, ct.ratio = 1,
          ia.adj = 1.0, alternative = c("two.sided", "one.sided"),
          tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{n}{Number of observations}
  \item{year}{Time point in which survival function is evaluated, default is 3}
  \item{s0}{Survival function at yaer = year under H0, default is 0.7}
  \item{hr}{Hazard ratio}
  \item{accrual}{Years of accrual, default is 3 years}
  \item{followup}{Years of followup, default is 3 years}
  \item{sig.level}{Significance level (Type I error probability)}
  \item{power}{Power of test (1 minus Type II error probability)}
  \item{ct.ratio}{Control / treatment ratio, default is 1}
  \item{ia.adj}{Ajdusted for interim analysis, default is 1 (no adjustment), ia.adj must > or = 1}
  \item{alternative}{Alternative, either two-sided or one-sided, default is two-sided}
  \item{tol}{tol for uniroot function}
}

\details{
These calculations use Poisson process and person-year under exponentioal distribution.

Exactly one of the parameters 'n','power' and 'hr' must be passed as NULL, and that parameter is determined from the others. 
Notice that the last one has non-NULL default so NULL must be explicitly passed if you want to compute it.
}
\value{
Object of class 'power.htest', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.

  \item{n.trt}{Sample size for the treatment arm}
  \item{n.ctl}{Sample size for control arm}
  \item{s1}{Survival rate for the treatment arm}
  \item{s0}{Survival rate for control arm}
  \item{delta}{Survival benefit at year = year}
  \item{n.annual}{Annual accrual rate}
  \item{n.event}{Total number of events required}
}
\references{
1. Lawless, J. Statistical Models and Methods for Lifetime Data, John Wiley and Sons, 1982.
2. Chow, S. C., Shao, J. and Wang, H. S. Sample szie calculations in clinical research (second edition, 2008). page 177. equation 7.3.6 
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
\code{\link{power.surv.test}}, 
\code{\link{pwr.2p.test}}
}
\examples{
##
pw = power.2surv.test(power = 0.8, s0 = 0.83, hr = 0.75, year = 5)
print(pw, digits = 4)
}
\keyword{Power} 
\keyword{Sample size}
\keyword{Two-sample}
\keyword{Survival}
\keyword{Clinical trial}
