\name{power.fdr.tt}
\alias{power.fdr.tt}

\title{
Power and sample size for False Discovery Rate (FDR) with t-tests
}
\description{
  Compute power and sample size of test or determine the effect size under alternative hypothesis to obtain target power for false discovery rate(same as power.anova.test).
}
\usage{
  power.fdr.tt(n1=NULL, cc.ratio=1, pi0=0.99, fdr=0.05, delta=NULL, 
              power = NULL, sigma=1, tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{n1}{Number of cases}
  \item{cc.ratio}{Control / case ratio, default is 1}
  \item{fdr}{False Discovery Rate (FDR)}
  \item{power}{Power of test (1 minus Type II error probability)}
  \item{pi0}{Percentage of null, default is 0.99}
  \item{delta}{Difference under althernative hypothesis}
  \item{sigma}{Variance of distribution}
  \item{tol}{tol for uniroot function}
}

\details{
Calculate power and sample size of study that requires control of False Discovery Rate (for example, 
analysis of RNA Micro-Array data.
The algorithm is based on methods proposed by Liu and Hwang (2007).
  
Exactly one of the parameters 'n1','power' and 'delta' must be passed as NULL, and that parameter is determined from the others. 
Notice that the last one has non-NULL default so NULL must be explicitly passed if you want to compute it.
}
\value{
  Object of class 'power.htest', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.
  \item{case.control}{A vector for number of cases and controls}
  \item{n.total}{Total number of sample size for the study}
  \item{pi0}{Percentage of null distributions}
  \item{sigma} {Varaince of the distribution}
  \item{fdr}{False Discovery Rate}
  \item{Power} {Power of the test}
  \item{delta} {Difference of means under althernative hypothesis}
  \item{c.reject} {Critical value to reject the null hypothesis}
  \
}
\references{
Liu, P. and Hwang J T.(2007). Quick calculation for sample size while controlling false
discovery rate with application to microarray analysis. Bioinformatics, 23, 739-746.
}
\author{
  Bingshu E. Chen (bingshu.chen@queensu.ca)
}
\note{
  'uniroot' is used to solve power equation for unknowns, so you may see errors from it, notably about inability to bracket the root when invalid arguments are given.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
  \seealso{
    \code{\link{power.pga.snp}}, 
    \code{\link{power.surv.test}},
    \code{\link{power.t.test}}, 
    \code{\link{pwr.2p.test}}
  }
\examples{
# Example from paper by Menash, Ronenberg and Chen (2008).
power.fdr.tt(n1 = 45, delta = 0.5, pi0 = 0.9, fdr = 0.1)
}
\keyword{ Power} 
\keyword{ Sample size}
\keyword{ False Discoverary Rate}
\keyword{ Micro-Array data analysis}
\keyword{ t-test}
