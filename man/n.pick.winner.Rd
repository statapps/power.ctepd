\name{n.pick.winner}
\alias{n.pick.winner}

\title{
Sample size calculation for Simon's pick winner design
}
\description{
Calculate sample size of Simon's pick winner design or determine probability
of alternative hypothesis to obtain target power.
}
\usage{
n.pick.winner(n = NULL, k = 2, p0=NULL, p1=NULL, power=NULL, 
                tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{n}{Number of sample size per arm}
  \item{k}{Number of arms, k > or = 2, default is k = 2}
  \item{p0}{Response rate of the lower efficacy drug}
  \item{p1}{Response rate of the best drug, p0 < p1}
  \item{power}{Probability of choosing the winner arm, power = 0.90 or higher is recommended}
  \item{tol}{tol for uniroot function}
}

\details{
These calculations use method proposed by R. Simon for pick winner design
for randomized phase II clinical trials.

Exactly one of the parameters 'n','power' and 'p1' must be passed as NULL, and that parameter is determined from the others.
Notice that the last one has non-NULL default so NULL must be explicitly passed if you want to compute it.

Since pick winner design is not a hypothesis test, there is not formal 
comparison between arms. 
Therefore, the type I error rate (sig.level) does not apply in this 
situation. 

'pselect' in package 'clinfun' also provide the probability of slection the winner arm but do not provide a direct solution for sample size.
}
\value{
Object of class 'power.htest', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.

  \item{n}{Number of sample size per arm}
  \item{arm}{Number of arms, the same as input k}
  \item{n.total}{Number of total sample size = n*arms}
}
\references{
Richard Simon, Robert E. Wittes and Susan S. Ellenberg. (1985).
Randomized phase II clinical Trials. Cancer Treat Rep. 
vol. 69, page 1375-1381.
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
\code{\link{pwr.p.test}}
}
\examples{
## Example: sample size for a three arms study
n.pick.winner(power = 0.95, p0 = 0.2, p1 = 0.35, k =3 )
}
\keyword{ Pick winner design }
\keyword{ Sample size}
\keyword{ Randomized phase II clinical trials}
