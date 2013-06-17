\name{OR}
\alias{OR}
\title{
OR(odd ratio) for risk_allele
}
\description{
OR(odd ratio) for risk_allele.
}
\usage{
OR(case_allele1, case_allele2, control_allele1, control_allele2)
}
\arguments{
  \item{case_allele1}{
a non-negative integer, the frequency of allele1 in case samples. 
}
  \item{case_allele2}{
a non-negative integer, the frequency of allele2 in case samples. 
}
  \item{control_allele1}{
a non-negative integer, the frequency of allele1 in control samples. 
}
  \item{control_allele2}{
a non-negative integer, the frequency of allele2 in control samples. 
}
}
\details{
The odds ratio is a measure of effect size, describing the strength of association or non-independence
between two binary data values.
}
\value{
  \item{risk_allele }{risk_allele(OR>1).}
  \item{OR }{the value of OR.}
  \item{CI }{a matrix with one row, the 0.95 CI(Confidence Interval) of OR: lower limit OR, upper limit OR.}
}
\references{
Edwards, A.W.F.(1963): The measure of association in a 2x2 table.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
   \code{\link{OR.TradPerm}},
	\code{\link{OR.MCPerm}},
   \code{\link{Armitage}},
	\code{\link{Armitage.TradPerm}},
   \code{\link{Armitage.MCPerm}},
	\code{\link{chisq.test}},
   \code{\link{chisq.TradPerm}},
	\code{\link{chisq.MCPerm}},
   \code{\link{fisher.test}},
	\code{\link{fisher.TradPerm}},
	\code{\link{fisher.MCPerm}},
	\code{\link{meta}},
	\code{\link{meta.TradPerm}},
	\code{\link{meta.MCPerm}},
	\code{\link{permuteGenotype}},
   \code{\link{rhyper}},
	\code{\link{permuteGenotypeCount}},
	\code{\link{genotypeStat}}
}
\examples{
# case_allele1=23
# case_allele2=34
# control_allele1=26
# control_allele2=29
# OR(case_allele1, case_allele2, control_allele1, control_allele2)
}
\keyword{ OR }
