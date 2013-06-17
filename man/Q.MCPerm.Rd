\name{Q.MCPerm}
\alias{Q.MCPerm}
\title{
Correct p.value for Heterogeneity statistics Q in meta analysis by MCPerm method.
}
\description{
Correct p.value for Heterogeneity statistics Q in meta analysis by MCPerm method.
}
\usage{
Q.MCPerm(case_11, case_12, case_22, control_11, control_12, control_22, 
    model = "allele", method = "MH", repeatNum = 1000)
}
\arguments{
  \item{case_11}{
a non-negative integer, the frequency of genotype "allele1/allele1" in case samples. 
}
  \item{case_12}{
a non-negative integer, the frequency of genotype "allele1/allele2" in case samples. 
}
  \item{case_22}{
a non-negative integer, the frequency of genotype "allele2/allele2" in case samples. 
}
  \item{control_11}{
a non-negative integer, the frequency of genotype "allele1/allele1" in control samples. 
}
  \item{control_12}{
a non-negative integer, the frequency of genotype "allele1/allele2" in control samples. 
}
  \item{control_22}{
a non-negative integer, the frequency of genotype "allele2/allele2" in control samples. 
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

MCPerm details see \code{\link{chisq.MCPerm}}.
}
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
\seealso{
  \code{\link{meta.MCPerm}},
  \code{\link{meta.TradPerm}},
  \code{\link{Q.TradPerm}},
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
# data(MetaGenotypeCount)
## delete the first line which is the names for columns.
# temp=MetaGenotypeCount[-1,,drop=FALSE]
# result=Q.MCPerm(case_11=as.numeric(temp[,14]),case_12=as.numeric(temp[,16]),
	 # case_22=as.numeric(temp[,18]),control_11=as.numeric(temp[,15]),
	 # control_12=as.numeric(temp[,17]),control_22=as.numeric(temp[,19]),
	 # model="allele",method="MH",repeatNum=100000)
# result
}
\keyword{ Q }
\keyword{ Q.MCPerm }
\keyword{ MCPerm }
