\name{Armitage.TradPerm}
\alias{Armitage.TradPerm}
\title{
A permutation method for Armitage's trend test in case/control association study
}
\description{
A permutation method for Armitage's trend test in case/control association study.
}
\usage{
Armitage.TradPerm(genotypeLine, affectionLine, fromCol, naString, sep, repeatNum = 1000)
}
\arguments{
  \item{genotypeLine}{
a matrix with one row containing information of specified snp: basic information(e.g. SNP 
ID number,chromo,position) and genotype of observed individuals. See below for details.
}
  \item{affectionLine}{
a matrix having the same dimension with parameter 'genotypeLine' contain the affection status
(case or control) of each individual and other information. The affection status must be in 
same columns with the genotypes in parameter 'genotypeLine'. See below for details.
}
  \item{fromCol}{
a positive integer, the start column of genotype data in parameter 'genotypeLine'.
}
  \item{naString}{
a character string for NA values of genotype.
}
  \item{sep}{
character separator used to divide genotype between alleles "Allele1<sep>Allele2".
}
  \item{repeatNum}{
an integer(default 1000) specifying the number of replicates for permutation test.
}
}
\details{
   The Cochran-Armitage test for trend, is used in categorical data analysis when the aim is to assess
for the presence of an association between a variable with two categories and a variable with k
categories. It modifies the Pearson chi-squared test to incorporate a suspected ordering in the effects
of the k categories of the second variable. The trend test is ofen used as a genotype-based test for 
case/control genetic association studies.

   'Armitage.TradPerm' uses permutation test for multiple Armitage's trend test correction. Seealso \code{\link{chisq.TradPerm}}.
   
   The basic information of sepcified snp for 'genotypeLine' and the other information of individuals
for 'affectionLine' must be located in the matrix of previous columns;and also can not be included,thus 
fromCol=1. 

   The genotypes of the specified snp, the stored alleles is considered to be ordered, i.e. "C/T" is 
unequivalent to "T/C".

   The affection status can be character string or numeric, but the symbol for control must be in advantageous 
alphabetical order (e.g. control=0,case=1).
}
\value{
   \item{pValue }{the p value for the test.}
   \item{obsStatistic }{the statistic of Armitage's trend test for the true data.}
   \item{obsP }{the p value for Armitage's trend test of the true data.}
   \item{permStatistic }{a vector with 'repeatNum' elements, the statistic of Armitage's trend test for the permutation data.}
   \item{permP }{a vector with 'repeatNum' elements, the p value for Armitage's trend test for the permutation data.}   
}
\references{
William S Noble(Nat Biotechnol.2009): How does multiple testing correction work?

Armitage, P(1955): Tests for Linear Trends in Proportions and Frequencies.

statgen.org(2007): A derivation for Armitage's trend test for the 2x3 genotype table.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
When input the parameter 'naString' and 'sep', please make sure correct. Otherwise the result will be wrong. 
}
\seealso{
   \code{\link{OR}},
	\code{\link{OR.TradPerm}},
	\code{\link{OR.MCPerm}},
   \code{\link{Armitage}},
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
## import example data
# data(genotypeData)
## get the first line: affection state for samples
# data1=genotypeData[1,,drop=FALSE]
## get the second line: genotype data for a sepcifed snp
# data2=genotypeData[2,,drop=FALSE]
# Armitage.TradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",repeatNum=10000)

## matrix
# genotypeLine=matrix(c("rs12","AA","TT","TA","AA","TT","AA","AA"),nrow=1)
# affectionLine=matrix(c("Affection",1,1,1,0,0,0,0),nrow=1)
# Armitage.TradPerm(genotypeLine,affectionLine,fromCol=2,naString="NN",sep="",repeatNum=5)

## connect file
# datafile=file("F:/data.txt","r")
## get the affection status for samples, read the first line from "data.txt"
# dataLine1=readLins(datafile,n=1)
# dataLine1=t(unlist(strsplit(dataLine1,sep="")))
## get the genotype line for samples, read the second line form "data.txt"
# dataLine2=readLines(datafile,n=1)
# dataLine2=t(unlist(strsplit(dataLine2,sep="")))
# Armitage.TradPerm(dataLine2,dataLine1,fromCol=2,naString="NN",sep="",repeatNum=5)
}
\keyword{ Armitage.TradPerm }
\keyword{ TradPerm }
