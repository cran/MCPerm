\name{OR.MCPerm}
\alias{OR.MCPerm}
\title{
A Monte Carlo permutation method for multiple OR(odd ratio) test correction in case/control association study
}
\description{
A Monte Carlo permutation method for multiple OR(odd ratio) test correction in case/control association study.
}
\usage{
OR.MCPerm(case_allele1, case_allele2, control_allele1, control_allele2, repeatNum = 1000)
}
\arguments{
  \item{case_allele1   }{
a non-negative integer, the frequency of allele1 in case samples. 
}
  \item{case_allele2   }{
a non-negative integer, the frequency of allele2 in case samples. 
}
  \item{control_allele1}{
a non-negative integer, the frequency of allele1 in control samples. 
}
  \item{control_allele2}{
a non-negative integer, the frequency of allele2 in control samples. 
}
  \item{repeatNum}{
an integer(default 1000) specifying the number of replicates used in the Monte Carlo permutation.
}
}
\details{
The odds ratio is a measure of effect size, describing the strength of association or non-independence
between two binary data values.

"OR.MCPerm" simulates permutation method to get p value for risk_allele, by generating random numbers 
for the hypergeometric distribution based on the genotype statistic of original data. Seealso \code{\link{chisq.MCPerm}}.
}
\value{
  \item{risk_allele }{risk_allele.}
  \item{pValue }{the p value for the test.}
  \item{obsOR }{the OR value for the true data.}
  \item{permOR }{a vector with 'repeatNum' elements, the OR value for the simulation data.}
}
\references{
William S Noble(Nat Biotechnol.2009): How does multiple testing correction work? 

Edwards, A.W.F.(1963): The measure of association in a 2x2 table.

Edgington. E.S.(1995): Randomization tests, 3rd ed.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
   \code{\link{OR}},
	\code{\link{OR.TradPerm}},
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
# OR.MCPerm(23,34,26,29,repeatNum=100000)
}
\keyword{ OR.MCPerm }
\keyword{ MCPerm }
