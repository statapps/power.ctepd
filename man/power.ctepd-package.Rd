\name{power.ctepd-package}
\alias{power.ctepd-package}
\alias{power.ctepd-doc}
\docType{package}
\title{
Biomarker Threshold Models
}
\description{
This package fits a threshold regression model for predictive and prognostic biomarker effect 
for binary data and survival data with an unknown cutoff point of biomarker.
}
\details{
"power.ctepd" is a R package for Biomarker Threshold Models.
Please use the following steps to install 'bhthm' package:

1. First, you need to install the 'devtools' package. You can skip this step if you have 'devtools' installed in your R. Invoke R and then type

  install.packages("devtools")

2. Load the devtools package.

  library(devtools)

3. Install "power.ctepd" package with R commond

  install_github("statapps/power.ctepd")


"power.ctepd" uses different statistical methods to identify cut-point 
(thershold parameter) for the biomarker in either generalized 
linear models or Cox proportional hazards model.
}
\author{Bingshu E. Chen

Maintainer: Bingshu E. Chen <bingshu.chen@queensu.ca>
}
\references{
Chen, B. E., Jiang, W. and Tu, D. (2014). A hierarchical Bayes model for biomarker subset effects in clinical trials. Computational Statistics and Data Analysis. vol 71, page 324-334.
}
\keyword{power.ctepd}
\seealso{
\code{power.t.test, }
\code{pwr.p.test, }
\code{pwr.2p.test}
}
\examples{
ppt = power.pois.test(n = 75, s0 = 0.83, s1 = 0.92, sig.level = 0.1)
} 
