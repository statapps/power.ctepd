\name{conditional.power}
\alias{conditional.power}
\alias{cpower.surv}
\alias{cpower.tt}

\title{
Conditional power for sequential tests.
}
\description{
  Compute conditional power for Sequential tests in clinical trials.
}
\usage{
  conditional.power(tk=0, theta=0.5, sig.level=0.025, power = 0.8, tm = 0.5)

# Conditional power for survival analysis
#
  cpower.surv(hrk = 1.0, HR = 0.7, sig.level = 0.025, power = 0.8, tm = 0.5) 
#
#
# Conditional power for t-test
#
  cpower.tt(muk = 0.0, mu = 0.7, sig.level = 0.025, power = 0.8, 
      tm = 0.5, sigma = 1.0, cc.ratio = 1) 
#
#
}
\arguments{
  \item{cc.ratio}{Control / case ratio, default is 1}
  \item{sig.level}{Significance level (Type I error probability)}
  \item{power}{Power of test (1 minus Type II error probability)}
  \item{HR}{Target hazards ratio for survival outcome}
  \item{mu}{Target difference for continuous outcome}
  \item{sigma}{Standard error for continuous outcome}
  \item{theta}{Difference under althernative hypothesis}
  \item{muk}{Difference in mean for the k^th interim analysis}
  \item{tk}{Test statistics Tk for the k^th interim analysis}
  \item{hrk}{Observed harzards ratio for the k^th interim analysis}
  \item{tm}{Time of the k^th test (percentage of information)}
}

\details{
Calculate conditional power for sequential tests in clinical trials to provide information on 
whether or not to futile a study.
The algorithm is based on methods by Jennison and Turnbull (2000) Chapter 10, page 207-220.
}
\value{
  Object of class 'power.htest', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.
  \item{case.control}{A vector for number of cases and controls}
  \item{n.total}{Total number of sample size for the study}
  \item{cpower.low}{Lower conditional power for one-sided futility analysis}
  \item{cpower.up }{Upper conditional power for one-sided futility analysis}
  \item{cpower}{Conditional power for two-sided test}
  \item{Power}{Power of the test}
  \item{I.max}{Maximum information at the final analysis}
  \item{zk}{Standardized test statistics Zk at the k^th interim analysis}
  \item{time}{Time of the k^th interim analysis}
}
\references{
Jennison C. and Turnbull B. W. (2000). Group Sequential methods applications to clinical trials.
}
\author{
  Bingshu E. Chen (bingshu.chen@queensu.ca)
}
\note{
  'uniroot' is used to solve power equation for unknowns, so you may see errors from it, notably about inability to bracket the root when invalid arguments are given.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
  \seealso{
    \code{\link{power.surv.test}},
    \code{\link{power.t.test}}, 
    \code{\link{pwr.2p.test}}
  }
\examples{
# Example
conditional.power(tk=0, theta=0.5, sig.level=0.05, power = 0.8, tm = 0.5)
}
\keyword{ Power} 
\keyword{ Sample size}
\keyword{ Conditional power}
\keyword{ Sequential test}
