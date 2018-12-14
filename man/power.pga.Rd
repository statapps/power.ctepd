\name{power.pga.snp}
\alias{power.pga.snp}
\alias{power.pga.hap}

\title{
Power and sample size for case-control genetics studies
}
\description{
  Compute power and sample size of test or determine the relative risks (RR) under alternative hypothesis to obtain target power (same as power.anova.test).
}
\usage{
  power.pga.snp(n1=NULL, cc.ratio=1, p0=0.01, maf=0.05, daf=0.05, Dp=NULL, 
              R2=NULL, rr=NULL, rr2=NULL, sig.level=0.05, power = NULL, 
              edf=1, genModel=c('co.dom', 'co.dom2', 'dom', 'rec'),
              tol = .Machine$double.eps^0.25)

  power.pga.hap(n1=NULL, cc.ratio=1, p0=0.01, dhf=0.05, n.hap=2, 
              rr=NULL, rr2=NULL, sig.level=0.05, power = NULL, 
              genModel=c('co.dom', 'co.dom2', 'dom', 'rec'),
              tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{n1}{Number of cases}
  \item{cc.ratio}{Control / case ratio, default is 1}
  \item{sig.level}{Significance level (Type I error probability)}
  \item{power}{Power of test (1 minus Type II error probability)}
  \item{Dp}{D', default is 1, only specify one of Dp and R2}
  \item{R2}{R-square, default is 1}
  \item{p0}{Disease prevalence, default is 0.01}
  \item{maf}{Minor allele frequency, default is 0.05}
  \item{daf}{Disease allele frequency, default is 0.05}
  \item{dhf}{Disease haplotype frequency, default is 0.05}
  \item{n.hap}{Number of haplotypes, default is 2}
  \item{rr}{Relative risk (RR)}
  \item{rr2}{Second relative risk for co-dominant model}
  \item{genModel}{'co.dom': co-dominant mdoel (1 df), 'co.dom2', co-dominant model (2df), 'dom': dominant model, 'rec': recessive model}
  \item{edf}{edf: Effective number of tests, default is 1}
  \item{tol}{tol for uniroot function}
}

\details{
Calculate power and sample size of case-control study for genetics analysis (either SNP or haplotype) using a co-dominant, dominant or recessive disease susceptibility model.
The algorithm is based on methods proposed by Menash, Ronenberg and Chen (2008).
  
  Exactly one of the parameters 'n1','power' and 'rr' must be passed as NULL, and that parameter is determined from the others. 
  Notice that the last one has non-NULL default so NULL must be explicitly passed if you want to compute it.
}
\value{
  Object of class 'power.htest', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.
  
  \item{case.control}{A vector for number of cases and controls}
  \item{D_max}{Maximum disequilibrium D_max = min(p1(1-pd), pd(1-p1)), output for power.pga.snp}
  \item{D}{Non-scaled disequilibrium with D = Dp*D.max, or D = r*D.max, where D.max is the maximum disequilibrium, with D.max = min(p1(1-pd), pd(1-p1)). This is an output for power.pga.snp}
  \item{n.total}{Total number of sample size for the study}
  \item{power}{Power of the study}
  \item{rr}{Releative risk with one copy of disease allele or disease haplotype}
  \item{rr2}{Relative risk with two copies of diseaes allele or diseas haplotype}
  \item{sig.adjust}{Ajusted significant level based on number of tests for SNP analysis}
  \item{sig.level} {Significant level}
  \item{Penetrances}{The penetrance rates for 0, 1 and 2 copies of disease SNP or haplotype}
}
\references{
Menashe, I., Rosenberg, P. S. and Chen, B. E., (2008). PGA: power calculator for case-control genetic association analysis. BioMed Central (BMC) Genetics. vol. 9, page 36.
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
# Example from paper by Menash, Ronenberg and Chen (2008).
power.pga.snp(power = 0.9, rr = 2.0, genModel = 'co.dom', p0 = 0.07, edf=200)
power.pga.hap(power = 0.9, rr = 1.7, genModel = 'co.dom', p0 = 0.07, n.hap=12)
}
\keyword{ Power} 
\keyword{ Sample size}
\keyword{ Case control study}
\keyword{ Genetics}
\keyword{ SNP}
\keyword{ Haplotype}
