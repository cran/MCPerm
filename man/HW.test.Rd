\name{HW.test}
\alias{HW.test}
\title{
Hardy-weinberg equilibrium test
}
\description{
Hardy-weinberg equilibrium test for control samples.
}
\usage{
HW.test(count_11, count_12, count_22)
}
\arguments{
  \item{count_11}{
a non-negative integer or vector, the frequency of genotype "allele1/allele1" in control samples. 
}
  \item{count_12}{
a non-negative integer or vector, the frequency of genotype "allele1/allele2" in control samples. 
}
  \item{count_22}{
a non-negative integer or vector, the frequency of genotype "allele2/allele2" in control samples. 
}
}
\details{
Hardy-weinberg equilibrium test states that allele and genotype frequencies in a population will remain constant
from generation to generation in the absence of other evolutionary influences.
}
\value{
  \item{X2 }{statistic(s) for the Hardy-weinberg equilibrium test.}
  \item{p.value }{p.value(s) for the Hardy-weinberg equilibrium test.}
}
\references{
Emigh, T.H. (1980 Biometrics 36(4): 627-642): A comparison of tests for Hardy-weinberg equilibrium.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
Hardy-weinberg equilibrium test for control samples.
}
\seealso{
   \code{\link{OR}},
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
# case_11=c(12,23,34,33)
# case_12=c(4,34,53,4)
# case_22=c(7,5,23,9)
# HW.test(case_11,case_12,case_22)
}
\keyword{ HW }
\keyword{ HW.test }