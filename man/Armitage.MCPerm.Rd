\name{Armitage.MCPerm}
\alias{Armitage.MCPerm}
\title{
A Monte Carlo permutation method for Armitage's trend test in case/control association study
}
\description{
A Monte Carlo permutation method for Armitage's trend test in case/control association study.
}
\usage{
Armitage.MCPerm(case_11, case_12, case_22, control_11, control_12, control_22, 
    repeatNum = 1000)
}
\arguments{
  \item{case_11}{
a non-negative integer, the frequency of genotype "allele1/allele1" in case samples. 
}
  \item{case_12}{
a non-negative integer, the frequency of genotype "allele1/allele2" in case samples. 
}
  \item{case_22}{
a non-negative integer, the frequency of genotype "allele2/allele2" in case samples.
}
  \item{control_11}{
a non-negative integer, the frequency of genotype "allele1/allele1" in control samples.
}
  \item{control_12}{
a non-negative integer, the frequency of genotype "allele1/allele2" in control samples.
}
  \item{control_22}{
a non-negative integer, the frequency of genotype "allele2/allele2" in control samples
}
  \item{repeatNum}{
an integer(default 1000) specifying the number of replicates used in the Monte Carlo permutation.
}
}
\details{
The Cochran-Armitage test for trend, is used in categorical data analysis when the aim is to assess
for the presence of an association between a variable with two categories and a variable with k
categories. It modifies the Pearson chi-squared test to incorporate a suspected ordering in the effects
of the k categories of the second variable. The trend test is ofen used as a genotype-based test for 
case/control genetic association studies.

"Armitage.MCPerm" simulates permutation method to correct p value, by twice generating random numbers for the 
hypergeometric distribution based on the genotype statistic of original data. Seealso \code{\link{chisq.MCPerm}}.
}
\value{
   \item{pValue }{the p value for the test.}
   \item{obsStatistic }{the statistic of Armitage's trend test for the true data.}
   \item{obsP }{the p value for Armitage's trend test of the true data.}
   \item{permStatistic }{a vector with 'repeatNum' elements, the statistic of Armitage's trend test for the simulation data.}
   \item{permP }{a vector with 'repeatNum' elements, the p value for Armitage's trend test for the simulation data.}   
}
\references{
William S Noble(Nat Biotechnol.2009): How does multiple testing correction work?

statgen.org(2007): A derivation for Armitage's trend test for the 2x3 genotype table.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
   \code{\link{OR}},
	\code{\link{OR.TradPerm}},
	\code{\link{OR.MCPerm}},
   \code{\link{Armitage}},
   \code{\link{Armitage.TradPerm}},
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
## Armitage.MCPerm(case_11=4,case_12=1,case_22=1,control_11=3,
  # control_12=5,control_22=7,repeatNum=10000)
}
\keyword{ Armitage.MCPerm }
\keyword{ MCPerm }
