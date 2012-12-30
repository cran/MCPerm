\name{tradPerm}
\alias{tradPerm}
\title{
A permutation test for multiple test correlation in case-control association study
}
\description{
A permutation test for multiple test correlation in case-control association study.
}
\usage{
tradPerm(genotypeLine, affectionLine, fromCol, naString, sep, testMethod = "chisq", repeatNum = 1000)
}
\arguments{
  \item{genotypeLine}{
a matrix or data.frame with one row containing information of specified snp(rs): basic information(e.g. 
SNP ID number,chromo,position) and genotype of observed individuals. See below for details.
}
  \item{affectionLine}{
a matrix or data.frame having the same dimension with parameter 'genotypeLine' contain the affection 
status(case or control) of each individual and other information. The affection status must be in same columns
with the genotypes in parameter GenotypeLine. See below for details.
}
  \item{fromCol}{
integer,the start column of genotype data in parameter 'genotypeLine'.
}
  \item{naString}{
a character string for NA values of genotype
}
  \item{sep}{
character separator used to divide genotype encoding both alleles---"Allele1<sep>Allele2".
}
  \item{testMethod}{
character string,includes "chisq","fisher","OR","Armitage". The method is called to compute a p-value for case-
control association study using permutation method.
}
  \item{repeatNum}{
an integer specifying the number of replicates for permutation.
}
}
\details{
For case-control association study for snps, the permutation test proceeds as follows:
1) Combine the observations from all the samples;
2) Shuffle them and rearrangements of the labels(case/control) on the observed data;
3) Record the genotype frequency of case and control samples, respectively;
4) Calculate the statistic of interest;
5) Repeat many times(at least 1000) to obtain the distribution of the statistic;
6) Determine how often the resampled statistic of interest is as extreme as the observed value of the same statistic. 
 
The basic information of sepcified snp for GenotypeLine and the other information of individuals for 
affectionLine must be located in the matrix or data.frame of previous columns; and can not be included,thus
fromCol=1.
The affection status can be character string(e.g. case,control) or numeric.
When count genotype frequency of the specified snp, the stored alleles is considered to be unordered, 
i.e. "C/T" is equivalent to "T/C".
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
    \item{permOR }{the OR-value for the simulation data.}
	
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
When input the parameter 'naString' and 'sep',please make sure correct.Otherwise the result will be wrong.
If genotypeLine has no three genotype and testMethod="Armitage", the return value will be NA.
}
\seealso{
  \code{\link{genotypeStat}},
  \code{\link{permuteData}}, 
  \code{\link{permuteMatrix}},   
  \code{\link{OR}},   
  \code{\link{Armitage}}, 
  \code{\link{MCPerm}},
  \code{\link{chisq.test}}, 
  \code{\link{fisher.test}}, 
  \code{\link{rhyper}}
}
\examples{
## import example data
data(exampleData)
data1=exampleData[1,,drop=FALSE]
data2=exampleData[2,,drop=FALSE]
ORresult=tradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",testMethod="OR",repeatNum=10)

## read data file
# setwd("F:/")
# fileIndex=file("DataFile.txt","r")
# readLine=readLines(fileIndex,n=1)
# data1=t(unlist(strsplit(readLine," ")))
# readLine=readLines(fileIndex,n=1)
# data2=t(unlist(strsplit(readLine," ")))
# result=tradPerm(data2,data1,fromCol=2,naString="?_?",sep="_")
# close(fileIndex)
}
\keyword{permutation }

