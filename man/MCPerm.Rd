\name{MCPerm}
\alias{MCPerm}
\title{
A Monte Carlo permutation method for multiple test correlation in case-control association study
}
\description{
A Monte Carlo permutation method for multiple test correlation in case-control association study.
}
\usage{
MCPerm(genotypeStat, testMethod = "chisq", repeatNum = 1000)
}
\arguments{
  \item{genotypeStat}{
a matrix of genotype frequency with two rows indicating case and control and three column for genotypes
(i.e. AA,AG,GG). See note. 
}
  \item{testMethod}{
character string,includes "chisq","fisher","OR","Armitage". The method is called to compute a p-value for case-
control association study using a Monte Carlo permutation method.
}
  \item{repeatNum}{
an integer specifying the number of replicates for the Monte Carlo permutation method.
}
}
\details{
  Permutation tests exist for any test statistic, regardless of whether or not its distribution is known.
Thus the permutation test is widely considered the gold standard for accurate multiple testing correction. 
For case-control association study for snps, the permutation test proceeds as follows:
1) Combine the observations from all the samples;
2) Shuffle them and rearrangements of the labels(case/control) on the observed data;
3) Record the genotype frequency of case and control samples, respectively;
4) Calculate the statistic of interest;
5) Repeat many times(at least 1000) to obtain the distribution of the statistic;
6) Determine how often the resampled statistic of interest is as extreme as the observed value of the same statistic. 
  Obviously, for multiple test correlation in case-control association study for millions of snp, the traditional
method---permutation test is very computationally impractical. Thus propose an accurate, rapid and efficient 
method for multiple testing correction in genome-wide association studies---MCPerm. 
  Method---MCPerm generates the genotype frequency for rearranged case and control data by twice generating 
random numbers for the hypergeometric distribution, based on the genotype statistic of original data, taking the 
place of the step 2) and step 3) of the traditional method. And the genotype frequency distribution generating by
MCPerm is almost the same with permutation test, this simplified method greatly improves the efficiency of the 
permutation test and is faster. MCPerm method can be the perfect alternative to permutation test.

If the matrix of genotypeStat has two columns with zeros and testMethod="OR", the return value will be NA. 
}
\value{
If the testMethod is "chisq",the returned value is a list with the following components:
    \item{pValue }{the p-value for the test.}
    \item{obsStatistic }{the value of the chi-squared test statistic for the true data.}
    \item{obsP }{the p-value of the chisq.test for the true data.}
    \item{permStatistic }{a matrix with one row, the value of the chi-squared test statistic for simulation data.}
    \item{permP }{a matrix with one row,the p-value of the chisq.test for simulation data.}
  
If the testMethod is "fisher",the returned value is a list with the following components:
    \item{pValue }{the p-value for the test.}
    \item{obsP }{the p-value of the fisher.test for the true data.}
    \item{permP }{a matrix with one row,the p-value of the fisher.test for simulation data.}
  
If the testMethod is "OR",the returned value is a list with the following components:
    \item{pValue }{the p-value for the test.}
    \item{obsOR }{the OR-value for the true data.}
    \item{permOR }{a mtrix with one row,the OR-value for the simulation data.}

If the testMethod is "Armitage",the returned value is a list with the following components:
    \item{pValue }{the p-value for the test.}
    \item{obsStatistic }{the statistic of Armitage test for trend for the true data.}
    \item{obsP }{the p-value for Armitage test for trend of the true data.}	
	\item{permStatistic }{a matrix of one row,the statistic of Armitage test for trend for the simulation data.}
	\item{permP }{a matrix of one row,the p-value for Armitage test for trend for the simulation data.}
}
\references{
Edgington. E.S. (1995) Randomization tests, 3rd ed. New York: Marcel-Dekker
}
\author{
Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
If the matrix of genotypeStat has two columns with zeros and testMethod="OR", the return value will be NA. 
If the matrix of genotypeStat has no three columns or no two rows and testMethod="Armitage", the return value will be NA.
}
\seealso{
  \code{\link{genotypeStat}},
  \code{\link{permuteData}}, 
  \code{\link{permuteMatrix}},   
  \code{\link{OR}},  
  \code{\link{Armitage}},   
  \code{\link{tradPerm}},
  \code{\link{chisq.test}}, 
  \code{\link{fisher.test}}, 
  \code{\link{rhyper}}
}
\examples{
genotypeCount=matrix(c(34,0,45,13,0,9),nrow=2,byrow=TRUE)
MCPerm(genotypeCount,"chisq",repeatNum=10)
MCPerm(genotypeCount,"fisher",repeatNum=10)
MCPerm(genotypeCount,"OR",repeatNum=10)
MCPerm(genotypeCount,"Armitage",repeatNum=10)
}
\keyword{MCPerm}
