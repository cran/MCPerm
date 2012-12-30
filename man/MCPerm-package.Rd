\name{MCPerm-package}
\alias{MCPerm-package}
\docType{package}
\title{
A Monte Carlo permutation method for multiple test correlation in case-control association study
}
\description{
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
}
\details{
\tabular{ll}{
Package: \tab MCPerm\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2012-12-30\cr
License: \tab GPL-2\cr
}
}
\author{
Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\references{
Edgington. E.S. (1995) Randomization tests, 3rd ed. New York: Marcel-Dekker
}
\keyword{MCPerm }
\seealso{
  \code{\link{genotypeStat}},
  \code{\link{permuteData}}, 
  \code{\link{permuteMatrix}},   
  \code{\link{OR}},  
  \code{\link{Armitage}},   
  \code{\link{tradPerm}},
  \code{\link{MCPerm}},
  \code{\link{chisq.test}}, 
  \code{\link{fisher.test}}, 
  \code{\link{rhyper}}
}
\examples{
## example1-----genotypeStat-----------
# import example data
data(exampleData)
data1=exampleData[1,,drop=FALSE]
data2=exampleData[2,,drop=FALSE]
result=genotypeStat(data2,data1,fromCol=2,naString="?_?",sep="_")

## example2-----permuteData-----------
dataLine=matrix(c("rs12","1","2","3","4","5"),nrow=1)
newData=permuteData(dataLine=dataLine,fromCol=2)

## example3-----permuteMatrix-----------
rawMatrix=matrix(c(24,7,8,34),nrow=2)
newMatrix=permuteMatrix(rawMatrix)

## example4-----OR-----------
ORvalue=OR(20,30,10,60)

## example5-----tradPerm-----------
data(exampleData)
data1=exampleData[1,,drop=FALSE]
data2=exampleData[2,,drop=FALSE]
ORresult=tradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",testMethod="OR",repeatNum=10)

## example6-----MCPerm-----------
genotypeCount=matrix(c(34,0,45,13,0,9),nrow=2,byrow=TRUE)
result=MCPerm(genotypeCount,"chisq",repeatNum=10)
}
