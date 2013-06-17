\name{OR.TradPerm}
\alias{OR.TradPerm}
\title{
A permutation test for multiple OR(odd ratio) test correction in case/control association study
}
\description{
A permutation test for multiple OR(odd ratio) test correction in case/control association study.
}
\usage{
OR.TradPerm(genotypeLine, affectionLine, fromCol, naString, sep, repeatNum = 1000)
}
\arguments{
  \item{genotypeLine}{
a matrix with one row containing information of specified snp: basic information(e.g. SNP ID number,
chromo,position) and genotype of observed individuals. See below for details.
}
  \item{affectionLine}{
a matrix having the same dimension with parameter 'genotypeLine' contain the affection status(case or 
control) of each individual and other information. The affection status must be in same columns with 
the genotypes in parameter 'genotypeLine'. See below for details.
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
The odds ratio is a measure of effect size, describing the strength of association or non-independence
between two binary data values.

For case/control association study for snps, the permutation test proceeds as follows: 

1) Combine the observations from all the samples; 

2) Shuffle them and rearrangements of the labels(case/control) on the observed data; 

3) Record the allele frequency of case and control samples, respectively;

4) Calculate the statistic of interest; 

5) Repeat many times(at least 1000) to obtain the distribution of the statistic; 

6) Determine how often the resampled statistic of interest is as extreme as the observed
value of the same statistic.

   The basic information of sepcified snp for 'genotypeLine' and the other information of individuals
for 'affectionLine' must be located in the matrix of previous columns; and can not be included,thus 
fromCol=1. 

   The genotypes of the specified snp, the stored alleles is considered to be ordered, i.e. "C/T" is 
unequivalent to "T/C".

   The affection status can be character string or numeric, but the symbol for control must be in advantageous 
alphabetical order (e.g. control=0,case=1).
}
\value{
  \item{risk_allele }{risk_allele(OR>1).}
  \item{pValue }{the p value for the test.}
  \item{obsOR }{the OR value for the true data.}
  \item{permOR }{a vector with 'repeatNum' elements, the OR value for the permutation data.}
}
\references{
William S Noble(Nat Biotechnol.2009): How does multiple testing correction work? 

Edwards, A.W.F.(1963): The measure of association in a 2x2 table.

Edgington. E.S.(1995): Randomization tests, 3rd ed.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
When input the parameter 'naString' and 'sep', please make sure correct. Otherwise the result will be wrong.
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
## import example data
# data(genotypeData)
## get the first line: affection state for samples
# data1=genotypeData[1,,drop=FALSE]
## get the second line: genotype data for a sepcifed snp
# data2=genotypeData[2,,drop=FALSE]
# OR.TradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",repeatNum=5)

## matrix
# genotypeLine=matrix(c("rs12","AA","TT","TA","AA","TT","AA","AA"),nrow=1)
# affectionLine=matrix(c("Affection",1,1,1,0,0,0,0),nrow=1)
# OR.TradPerm(genotypeLine,affectionLine,fromCol=2,naString="NN",sep="",repeatNum=5)

## connect file
# datafile=file("F:/data.txt","r")
## get the affection status for samples, read the first line from "data.txt"
# dataLine1=readLins(datafile,n=1)
# dataLine1=t(unlist(strsplit(dataLine1,sep="")))
## get the genotype line for samples, read the second line form "data.txt"
# dataLine2=readLines(datafile,n=1)
# dataLine2=t(unlist(strsplit(dataLine2,sep="")))
# OR.TradPerm(dataLine2,dataLine1,fromCol=2,naString="NN",sep="",repeatNum=1000)
}
\keyword{ OR.TradPerm }
\keyword{ TradPerm }
