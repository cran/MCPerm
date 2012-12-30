\name{genotypeStat}
\alias{genotypeStat}
\title{
Statistical Allele and Genotype Frequency of the specified snp(rs#)
}
\description{
Statistical allele and genotype frequency of the specified snp(rs#) for case-control data.
}
\usage{
genotypeStat(genotypeLine, affectionLine, fromCol, naString, sep)
}
\arguments{
  \item{genotypeLine}{
a matrix or data.frame with one row containing information of specified snp(rs#): basic information(e.g. 
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
a character string for NA values of genotype.
}
  \item{sep}{
character separator used to divide genotype encoding both alleles---"Allele1<sep>Allele2".
}
}
\details{
The basic information of sepcified snp for GenotypeLine and the other information of individuals for 
affectionLine must be located in the matrix or data.frame of previous columns; and can not be included,thus
fromCol=1.
The affection status can be character string(e.g. case,control) or numeric.
When count genotype frequency of the specified snp, the stored alleles is considered to be unordered, 
i.e. "C/T" is equivalent to "T/C".
}
\value{
The returned value of genotypeStat is a list with the following components:
  \item{alleleStat }{
A matrix: the number of rows equal to the number of affection status; the number of columns equal to the 
number of the alleles of specified snp. The row names of the matrix are symbols of the affection status
in parameter 'affectionLine'; the column names of the matrix are symbols of the alleles in parameter 
'genotypeLine'. 
}
  \item{genotypeStat }{
A matrix: the number of rows equal to the number of affection status; the number of columns equal to the 
number of the genotypes of specified snp. The row names of the matrix are symbols of the affection status
in parameter 'affectionLine'; the column names of the matrix are symbols of the genotypes in parameter 
'genotypeLine'. The column order in matrix is ascending order by symbols of the genotypes, and the last 
column is for NA values. See note.
}
}
\references{
Edgington. E.S. (1995) Randomization tests, 3rd ed. New York: Marcel-Dekker
}
\author{
Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
the parameter 'naString' must be right, otherwise the columns of the value genotypeStat will not be orderd in 
the described order in value genotypeStat that the last column will be not for NA values.
}
\seealso{
  \code{\link{permuteData}}, 
  \code{\link{permuteMatrix}},   
  \code{\link{OR}},   
  \code{\link{Armitage}},    
  \code{\link{tradPerm}},
  \code{\link{MCPerm}}
}
\examples{
## import example data
data(exampleData)
data1=exampleData[1,,drop=FALSE]
data2=exampleData[2,,drop=FALSE]
result=genotypeStat(data2,data1,fromCol=2,naString="?_?",sep="_")

## read data file
# setwd("F:/")
# fileIndex=file("DataFile.txt","r")
# readLine=readLines(fileIndex,n=1)
# data1=t(unlist(strsplit(readLine," ")))
# readLine=readLines(fileIndex,n=1)
# data2=t(unlist(strsplit(readLine," ")))
# result=genotypeStat(data2,data1,fromCol=2,naString="?_?",sep="_")
# close(fileIndex)
}
\keyword{genotypeStat}
\keyword{genotype}
\keyword{allele}