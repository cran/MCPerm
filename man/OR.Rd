\name{OR}
\alias{OR}
\title{
OR(odd ratio) for the disease-allele 
}
\description{
OR(odd ratio) for the disease-allele.
}
\usage{
OR(caseAllele1, caseAllele2, controlAllele1, controlAllele2)
}
\arguments{
  \item{caseAllele1}{
integer,the allele1 frequency in case samples.
}
  \item{caseAllele2}{
integer,the allele2 frequency in case samples.
}
  \item{controlAllele1}{
integer,the allele1 frequency in control samples.
}
  \item{controlAllele2}{
integer,the allele2 frequency in control samples.
}
}
\details{
The odds ratio is a measure of effect size, describing the strength of association or non-independence between 
two binary data values.
}
\value{
  \item{OR }{the value of OR.}
  \item{CI }{a matrix with one row,the 0.95 CI(Confidence Interval) of OR:lower limit OR,upper limit OR.}
}
\references{
Edwards, A.W.F. (1963). "The measure of association in a 2x2 table".
}
\author{
Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
each parameter cannot be zero.
}
\seealso{
  \code{\link{genotypeStat}},
  \code{\link{permuteData}}, 
  \code{\link{permuteMatrix}},  
  \code{\link{Armitage}},    
  \code{\link{tradPerm}},
  \code{\link{MCPerm}}
}
\examples{
ORvalue=OR(20,30,10,60)
}
\keyword{OR}
