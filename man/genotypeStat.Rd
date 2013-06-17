\name{genotypeStat}
\alias{genotypeStat}
\title{
Statistical Allele and Genotype Frequency of the specified snp
}
\description{
Statistical Allele and Genotype Frequency of the specified snp.
}
\usage{
genotypeStat(genotypeLine, affectionLine, fromCol, naString, sep)
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
a character string for NA value of genotype.
}
  \item{sep}{
character separator used to divide genotype between alleles "Allele1<sep>Allele2".
}
}
\details{
   The basic information of sepcified snp for 'genotypeLine' and the other information of individuals
for 'affectionLine' must be located in the matrix of previous columns; and also can not be included,thus 
fromCol=1. 

   The genotypes of the specified snp, the stored alleles is considered to be ordered, i.e. "C/T" is 
unequivalent to "T/C".

   The affection status can be character string or numeric, but the symbol for control must be in advantageous 
alphabetical order (e.g. control=0,case=1).
}
\value{
  \item{alleleCount }{
A vector with four elements: 'case_allele1','case_allele2','control_allele1',
'control_allele2'.
}
  \item{genotypeCount }{
A vector with eight elements: 'case_11','case_12','case_22','case_NA','control_11',
'control_12','control_22','control_NA'.
}
}
\references{
William S Noble(Nat Biotechnol.2009): How does multiple testing correction work?

Edgington. E.S.(1995): Randomization tests, 3rd ed.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
When input the parameter 'naString' and 'sep', please make sure correct. Otherwise the result will be wrong. 
}
\seealso{
   \code{\link{OR.TradPerm}},
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
	\code{\link{permuteGenotypeCount}}
}
\examples{
## import example data(data.frame)
# data(genotypeData)
## get the first line: affection state for samples
#data1=genotypeData[1,,drop=FALSE]
## get the second line: genotype data for a sepcifed snp
# data2=genotypeData[2,,drop=FALSE]

## Statistical allele and genotype frequency of the specified snp(rs#) for case-control data.
# result2=genotypeStat(data2,data1,fromCol=2,naString="?_?",sep="_")
## allele frequency for case and control samples
# alleleResult=result2$alleleCount
## genotype frequency for case and control samples
# genotypeReslut=result2$genotypeCount

## matrix
# genotypeLine=matrix(c("rs12","AA","AT","NN","AT","AT","AT"),nrow=1)
# affectionLine=matrix(c("Affection",0,0,1,1,0,0),nrow=1)
# fromCol=2
# naString="NN"
# sep=""
# genotypeStat(genotypeLine, affectionLine, fromCol, naString, sep)

## connect file
# datafile=file("F:/data.txt","r")
## get the affection status of samples: the first line from "data.txt"
# dataLine1=readLines(datafile,n=1)
# dataLine1=t(unlist(strsplit(dataLine1,sep="")))
## get the genotypes of samples: the second line from "data.txt"
# dataLine2=readLines(datafile,n=1)
# dataLine2=t(unlist(strsplit(dataLine2,sep="")))
# genotypeStat(dataLine2,dataLine1,fromCol=2,naString="NN",sep="")
}
\keyword{ genotypeStat }
