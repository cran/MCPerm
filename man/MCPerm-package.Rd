\name{MCPerm-package}
\alias{MCPerm-package}
\alias{MCPerm}
\docType{package}
\title{
A Monte Carlo permutation method for multiple test correction in case/control association study
}
\description{
Permutation tests exist for any test statistic, regardless of whether or not its distribution is known.
Thus the permutation test is widely considered the gold standard for accurate multiple testing correction. 

    For case/control association study for snps, the permutation test proceeds as follows: 
	
  1) Combine the observations from all the samples; 
  
  2) Shuffle them and rearrangements of the labels(case/control) on the observed data;
  
  3) Record the genotype frequency of case and control samples, respectively; 
  
  4) Calculate the statistic of interest; 
  
  5) Repeat many times(at least 1000) to obtain the distribution of the statistic; 
  
  6) Determine how often the resampled statistic of interest is as extreme as the observed value of the 
same statistic. 

    Obviously, for multiple test correction in case/control association study for millions of snp, the 
traditional method---permutation test is very computationally impractical. Thus propose an accurate, rapid 
and efficient method for multiple testing correction in genome-wide association studies---MCPerm. 

    Method---MCPerm generates the genotype frequency for rearranged case and control data by twice generating random 
numbers for the hypergeometric distribution, based on the genotype statistic of original data, taking the place
of the step 2) and step 3) of the traditional method. And the genotype frequency distribution generating by 
MCPerm is almost the same with permutation test, this simplified method greatly improves the efficiency of 
the permutation test and is faster. MCPerm method can be the perfect alternative to permutation test.
}
\details{
\tabular{ll}{
Package: \tab MCPerm\cr
Type: \tab Package\cr
Version: \tab 1.1.3\cr
Date: \tab 2013-01-27\cr
License: \tab GPL-2\cr
}
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\references{
William S Noble(Nat Biotechnol.2009): How does mutiple testing correction work?

Edgington. E.S.(1995): Randomization tests, 3rd ed.
}
\keyword{ MCPerm }
\seealso{
    \code{\link{Armitage}},
	\code{\link{Armitage.TradPerm}},
    \code{\link{Armitage.MCPerm}},
    \code{\link{OR}},
    \code{\link{OR.TradPerm}},
	\code{\link{OR.MCPerm}},
	\code{\link{permuteData}},
	\code{\link{permuteGenotype}},
	\code{\link{genotypeStat}},
    \code{\link{chisq.TradPerm}},
	\code{\link{chisq.MCPerm}},
	\code{\link{fisher.TradPerm}},
	\code{\link{fisher.MCPerm}},
	\code{\link{rhyper}},
	\code{\link{chisq.test}},
	\code{\link{fisher.test}}
}
\examples{
## example1-----genotypeStat-----------
# import example data
data(exampleData)
# get the first line: affection state for samples
data1=exampleData[1,,drop=FALSE]
# get the second line: genotype data for a sepcifed snp
data2=exampleData[2,,drop=FALSE]
# Statistical allele and genotype frequency of the specified snp for case-control data.
result=genotypeStat(data2,data1,fromCol=2,naString="?_?",sep="_")
genotypeCount=result$genotypeCount
alleleCount=result$alleleCount

## example2-----permuteData-----------
# a matrix with 1 row
dataLine=matrix(c("rs12","1","2","3","4","5"),nrow=1)
# permutate the elements of the matrix
newData=permuteData(dataLine=dataLine,fromCol=2)

## example3-----permuteGenotype-----------
newMatrix=permuteGenotype(case_11=24,case_12=34,case_22=45,control_11=23,control_12=45,control_22=34,n=5)

## example4-----OR-----------
# OR(odd ratio) for the risk-allele
ORvalue=OR(case_allele1=20,case_allele2=30,control_allele1=10,control_allele2=60)

## example5-----OR.TradPerm--------
# import example data
data(exampleData)
# get the first line: affection state for samples
data1=exampleData[1,,drop=FALSE]
# get the second line: genotype data for a sepcifed snp
data2=exampleData[2,,drop=FALSE]
result1=OR.TradPerm(genotypeLine=data2,affectionLine=data1,fromCol=2,naString="?_?",sep="_",repeatNum=5)
risk_allele=result1$risk_allele
p=result1$pValue
obsOR=result$OR

## example6-----OR.MCPerm----------
OR.MCPerm(case_allele1=34,case_allele2=23,control_allele1=27,control_allele2=45,repeatNum=5)

## example7-----Armitage-------
Armitage(case_11=23,case_12=45,case_22=12,control_11=27,control_12=12,control_22=45)

## example8-----Armitage.TradPerm-----
# import example data
data(exampleData)
# get the first line: affection state for samples
data1=exampleData[1,,drop=FALSE]
# get the second line: genotype data for a sepcifed snp
data2=exampleData[2,,drop=FALSE]
Armitage.TradPerm(genotypeLine=data2,affectionLine=data1,fromCol=2,naString="?_?",sep="_",repeatNum=5)

## example9----Armitage.MCPerm--------
Armitage.MCPerm(case_11=23,case_12=45,case_22=12,control_11=27,control_12=12,control_22=45,repeatNum=5)

## example10---chisq.TradPerm---------
# import example data
data(exampleData)
# get the first line: affection state for samples
data1=exampleData[1,,drop=FALSE]
# get the second line: genotype data for a sepcifed snp
data2=exampleData[2,,drop=FALSE]
chisq.TradPerm(genotypeLine=data2,affectionLine=data1,fromCol=2,naString="?_?",sep="_",repeatNum=5)

## example11---chisq.MCPerm--------
case_11=23
case_12=45
case_22=12
control_11=27
control_12=12
control_22=45
chisq.MCPerm(23,45,12,27,12,45,repeatNum=5)

## example12---fisher.TradPerm------
# import example data
data(exampleData)
# get the first line: affection state for samples
data1=exampleData[1,,drop=FALSE]
# get the second line: genotype data for a sepcifed snp
data2=exampleData[2,,drop=FALSE]
fisher.TradPerm(genotypeLine=data2,affectionLine=data1,fromCol=2,naString="?_?",sep="_",repeatNum=5)

## example13---fisher.MCPerm-------
fisher.MCPerm(23,45,12,27,12,45,repeatNum=5)
}
