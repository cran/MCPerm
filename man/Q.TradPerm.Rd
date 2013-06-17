\name{Q.TradPerm}
\alias{Q.TradPerm}
\title{
Correct p.value for Heterogeneity statistics Q in meta analysis by TradPerm method
}
\description{
Correct p.value for Heterogeneity statistics Q in meta analysis by TradPerm method.
}
\usage{
Q.TradPerm(genotypeData, affectionData, split, sep, naString, model = "allele", 
    method = "MH", repeatNum = 1000)
}
\arguments{
  \item{genotypeData}{
a matrix with one column and multiple rows, each row contains genotype data for case and control 
samples of certain study. Note the field separtor of each line must be same, and same with parameter
'affectionData'.
}
  \item{affectionData}{
a matrix with one column and multiple rows, each row contains the affection stats of case and control samples of 
certain study which must correspond to 'genotypeData'. Note the field separtor of each line must be same,and same 
with parameter 'genotypeData'.
}
  \item{split}{
the field separator character, which separates elements on each line of the parameter 'genotypeData' and 
'affectionData'. 'Split' and 'sep' cannot be same.
}
  \item{sep}{
character separator used to divide genotype between alleles "Allele1<sep>Allele2" in parameter 'genotypeData'.
'Split' and 'sep' cannot be same.
}
  \item{naString}{
a character string for NA values of genotype data in parameter 'genotypeData'.
}
  \item{model}{
a character string indicating the type of model("allele","dominant" or "recessive") supplied to the data.
The risk allele(see details) is marked as allele1. The allele model indicates allele1 versus allele2, the 
dominant model indicates <allele1/allele1> + <allele1/allele2> versus <allele2/allele2>,
the recessive model indicates <allele1/allele1> versus <allele1/allele2> + <allele2/allele2>.
Default value is 'allele'.
}
  \item{method}{
a character string indicating the method('Inverse','MH' or 'Peto') to calculate Q value. Default value is 
'MH'. See details.
}
  \item{repeatNum}{
an integer specifying the number of replicates used in the Monte Carlo permutation. Default value is 1000.
}
}
\details{
Allele 1 and allele 2 to each study have OR values. The risk allele is the allele which 
the number of studies which OR>1 more than half of the number of all studies.

Q value fit chi square distribution,the number of studies in meta analysis minus one is
degree of freedom of Q value.

TradPerm details see \code{\link{chisq.TradPerm}}.}
\value{
  \item{risk_allele }{the symbol of risk allele.See details.}
  \item{Q }{the Q statistics for observed meta data.}
  \item{Qp }{the p value for Q,See details.}
  \item{corrected_Qp }{the percentage of less than Qp value in simulative Qp values.}
}
\references{
Julian P.T.Higgins, Simon G.Thompson(Statistics in Medicine,2002): Quantifying heterogeneity in a meta-analysis.

Julian P.T.Higgins, Simon G.Thompson, Jonathan J Deeks(BMJ,2003):Measuring inconsistency in meta-analyses.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
'Split':the field separtor of each line for parameter 'genotypeData' and 'affectionData' must be same.
'Split' and 'sep' cannot be same.
}
\seealso{
  \code{\link{meta.MCPerm}},
  \code{\link{meta.TradPerm}},
  \code{\link{Q.MCPerm}},
  \code{\link{I2.TradPerm}},
  \code{\link{I2.MCPerm}},
  \code{\link{chisq.MCPerm}},
  \code{\link{chisq.TradPerm}},
  \code{\link{VS.Hist}},
  \code{\link{VS.KS}},
  \code{\link{VS.Allele.Hist}},
  \code{\link{VS.Genotype.Hist}},
  \code{\link{PermMeta.LnOR.Hist}},
  \code{\link{PermMeta.Hist}}
}
\examples{
## import data
# data(MetaGenotypeData)
## delete first line which contains the names of each column
# temp=MetaGenotypeData[-1,];
# rowNum=nrow(temp)
# gen=matrix(0,nrow=rowNum,ncol=1);
# aff=matrix(0,nrow=rowNum,ncol=1);
# for(j in 1:rowNum){
	 # gen[j,]=paste(temp[j,14],temp[j,15],sep=" ");
	 # case_num=length(unlist(strsplit(temp[j,14],split=" ")));
	 # control_num=length(unlist(strsplit(temp[j,15],split=" ")));
	 # case_aff=paste(rep(2,case_num),collapse=" ");
	 # control_aff=paste(rep(1,control_num),collapse=" ");
	 # aff[j,]=paste(case_aff,control_aff,sep=" ");
# }
# result=Q.TradPerm(gen,aff,split=" ",sep="/",naString="-",
    # model="allele",method="MH",repeatNum=1000) 
# result
}
\keyword{ Q }
\keyword{ Q.TradPerm }
\keyword{ TradPerm }
