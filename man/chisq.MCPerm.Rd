\name{chisq.MCPerm}
\alias{chisq.MCPerm}
\title{
A Monte Carlo permutation method for multiple chisq.test correction in case/control association study
}
\description{
A Monte Carlo permutation method for multiple chisq.test correction in case/control association study.
}
\usage{
chisq.MCPerm(case_11, case_12, case_22, 
    control_11, control_12, control_22, repeatNum = 1000)
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
a non-negative integer, the frequency of genotype "allele2/allele2" in control samples. 
}
  \item{repeatNum}{
an integer(default 1000) specifying the number of replicates used in the Monte Carlo permutation.
}
}
\details{
Permutation tests exist for any test statistic, regardless of whether or not its distribution is known.
Thus the permutation test is widely considered the gold standard for accurate multiple testing correction. 

    For case/control association study for snps, the permutation test proceeds as follows: 
	
  1) Combine the observations from all the samples; 
  
  2) Shuffle them and rearrangements of the labels(case/control) on the observed data; 
  
  3) Record the genotype frequency of case and control samples, respectively; 
  
  4) Calculate the statistic of interest; 
  
  5) Repeat many times(at least 1000) to obtain the distribution of the statistic; 
  
  6) Determine how often the resampled statistic of interest is as extreme as the observed value of the same 
statistic. 

    Obviously, for multiple test correction in case/control association study for millions of snp, the 
traditional method---permutation test is very computationally impractical. Thus propose an accurate, rapid 
and efficient method for multiple testing correction in genome-wide association studies---MCPerm. 

   Method---MCPerm generates the genotype frequency for rearranged case and control data by twice generating 
random numbers for the hypergeometric distribution, based on the genotype statistic of original data, taking 
the place of the step 2) and step 3) of the traditional method. And the genotype frequency distribution generating 
by MCPerm is almost the same with permutation test, this simplified method greatly improves the efficiency of 
the permutation test and is faster. MCPerm method can be the perfect alternative to permutation test.
}
\value{
  \item{pValue }{the p value for the test.}
  \item{obsStatistic }{the value of the chi-squared test statistic for the true data.}
  \item{obsP }{the p value of the chisq.test for the true data.}
  \item{permStatistic }{a matrix with one row and 'repeatNum' columns, the value of the chi-squared test statistic for simulation data.}
  \item{permP }{a matrix with one row and 'repeatNum' columns, the p value of the chisq.test for simulation data.}  
}
\references{
William S Noble(Nat Biotechnol.2009): How does multiple testing correction work?

Edgington. E.S.(1995): Randomization tests, 3rd ed.
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
   \code{\link{Armitage.MCPerm}},
	\code{\link{chisq.test}},
	\code{\link{chisq.TradPerm}},
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
# case_11=34
# case_12=0
# case_22=16
# control_11=14
# control_12=0
# control_22=13
# chisq.MCPerm(case_11,case_12,case_22,control_11,control_12,control_22,repeatNum=10000)
}
\keyword{ chisq.MCPerm }
\keyword{ MCPerm }
