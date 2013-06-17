\name{genotypeData}
\alias{genotypeData}
\docType{data}
\title{
Genotype Data from GWA16
}
\description{
Genotype Data from GWA16.
}
\usage{data(genotypeData)}
\format{
  A matrix with 10 snps with 2062 samples.
  The format is:
  first row$: affection state for 2062 samples;
  other rows$: genotype data for each snp, e.g. A_A A_G G_G ?_? ...
}
\source{
GAW16 Data from the North American Rheumatoid Arthritis Consortium (NARAC)
}
\references{
MacCluer JW, Cupples LA and Almasy L (eds) Genetic Analysis Workshop 16: Approaches to Analysis of Genome-Wide Data. Genetic Epidemiology 33 (Suppl 1), S1-S110 (2009). 
}
\examples{
## import example data(data.frame)
# data(genotypeData)
## get the first line: affection state for samples
# data1=genotypeData[1,,drop=FALSE]
## get the second line: genotype data for a sepcifed snp
# data2=genotypeData[2,,drop=FALSE]

## Statistical allele and genotype frequency of the specified snp(rs#) for case-control data.
# result2=genotypeStat(data2,data1,fromCol=2,naString="?_?",sep="_")
## allele frequency for case and control samples
# alleleResult=result2$alleleCount
## genotype frequency for case and control samples
# genotypeReslut=result2$genotypeCount

## use permutation test to do multiple test correction in case/control association study
## return the correlated p_value and other information, see 'chisq.TradPerm'.
# result1=chisq.TradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",repeatNum=1000)
# p1=result1$pValue
# result2=fisher.TradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",repeatNum=1000)
# p2=result2$pValue
# result3=Armitage.TradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",repeatNum=1000)
# p3=result3$pValue
# result4=OR.TradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",repeatNum=1000)
# risk_allele=result4$risk_allele
# p4=result4$pValue
}
\keyword{GWA16}
