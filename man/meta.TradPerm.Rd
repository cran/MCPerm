\name{meta.TradPerm}
\alias{meta.TradPerm}
\title{
Meta analysis corrected by permutation test
}
\description{
Meta analysis corrected by permutation test.
}
\usage{
meta.TradPerm(genotypeData, affectionData, split, sep, naString, 
    model = "allele", fixed_method = "MH", random_method = "DL", 
	 Qp_alpha = 0.01, repeatNum = 1000)
}
\arguments{
  \item{genotypeData}{
a matrix with one column and multiple rows, each row contains genotype data for case and control 
samples of certain study. Note the field separtor of each line must be same, and same with parameter
'affectionData'.}
  \item{affectionData}{
a matrix with one column and multiple rows, each row contains the affection stats of case and control samples of 
certain study which must correspond to 'genotypeData'. Note the field separtor of each line must be same,and same 
with parameter 'genotypeData'.
}
  \item{split}{
the field separator character, which separates elements on each line of the parameter 'genotypeData' and 'affectionData'.
'Split' and 'sep' cannot be same.
}
  \item{sep}{
character separator used to divide genotype between alleles "Allele1<sep>Allele2" in parameter 'genotypeData'.
'Split' and 'sep' cannot be same.
}
  \item{naString}{
a character string for NA values of genotype data in parameter 'genotypeData'.
}
  \item{model}{
a character string indicating the type of model("allele"(default),"dominant" or "recessive") supplied to the data.
The risk allele(see details) is marked as allele1. The allele model indicates allele1 versus allele2, the 
dominant model indicates <allele1/allele1> + <allele1/allele2> versus <allele2/allele2>,
the recessive model indicates <allele1/allele1> versus <allele1/allele2> + <allele2/allele2>.
}
  \item{fixed_method}{
a character string indicating the method('Inverse','MH'(default) or 'Peto') to fit fixed-effects model.
}
  \item{random_method}{
a character string indicating the method to fit random-effects model: "HE", "DL"(default), "SJ", "ML", "REML", or "EB".
}
  \item{Qp_alpha}{
the threshold value(0~1)(default 0.01) to refuse the null hypothesis that all studies were evaluating the same effect.
Not rejecting the above hypothesis usually leads a meta-analysis to adopt a fixed-effects model.
The fixed-effects model assumes that the estimated effect sizes only differ by the sampling error. 
In contrast, if a significant Q-statistic P<Qp_alpha indicates heterogeneity across studies, a 
random-effects model should be adopted.
}
  \item{repeatNum}{
an integer(default 1000) specifying the number of replicates used in the Monte Carlo permutation.
}
}
\details{
Allele 1 and allele 2 to each study have OR values. The risk allele is the allele which 
the number of studies which OR>1 more than half of the number of all studies.

TradPerm details see \code{\link{chisq.TradPerm}}.
}
\value{
  \item{corrected_result }{matrix with 3 rows and 3 columns. The first column is Cochran's Q-statistics value; the second
     column is I2(=max(Q-df/Q,0), the degree of inconsistency across studies); the third column is the merged log odds ratio.
     The first row is the value for the 3 columns calculated by true data of the studies; the second row is the p value calculated 
     by true data of the studies; the third row is the p value calculated by permutation data, namely equal to the proportion of the  
     statistic of interest as extreme as the observed value of the same statistic. Note: the p value for I2 donot 
     calculate throngh meta-analysis, so the value is NA.  }
  \item{risk_allele }{the symbol of risk allele.See details.}
  \item{true_merged_LnOR }{the merged log odd ratio by meta analysis using observation data. }
  \item{true_merged_LnOR_VAR }{variation of the merged log odd ratio by meta analysis using observation data. }
  \item{true_merged_LnOR_p }{p value for the merged log odd ratio by meta analysis using observation data. }
  \item{true_merged_LnOR_ci.lb }{lower bound of the confidence intervals for the merged log odds ratio. }
  \item{true_merged_LnOR_ci.ub }{upper bound of the confidence intervals for the merged log odds ratio. }
  \item{study_num }{the number of studies in the meta analysis.}
  \item{sample }{a vector with 'study_num' elements, the sample size of each study.}
  \item{true_LnOR }{a vector of log odd ratio calculated by each study. }
  \item{true_VARLnOR }{a vector, the varition of log odd ratio calculated by each study. }
  \item{perm_case_11,perm_case_12,perm_case_22,perm_control_11,perm_control_12,perm_control_22 }{a matrix with 
    one row and 'repeatNum' columns, the freqency for genotype got by permutating data by MCPerm method. MCPerm
	 method details see \code{\link{chisq.MCPerm}}.}
  \item{perm_LnOR }{matirx with rows indicating studies and 'repeatNum' columns indicating log odds ratio for each permutation data.}
  \item{perm_VARLnOR }{matrix with rows indicating studies and 'repeatNum' columns indicating variance for the log odds ratio for each permuation data.}
  \item{perm_Qp }{a matrix with one row and 'repeatNum' columns, statistic Q value of heterogeneity for each permutation data.}
  \item{perm_I2 }{a matrix with one row and 'repeatNum' columns, statistic I2 value of heterogeneity for each permutation data.}
  \item{perm_merged_LnOR }{a matrix with one row and 'repeatNum' columns, merged log odd ratio for permutation data.}
  \item{perm_merged_VARLnOR }{a matrix with one row and 'repeatNum' columns, variation of merged log odd ratio for permuation data.}
  \item{perm_p }{a matrix with one row and 'repeatNum' columns, p value for merged log odd ratio of permuation data.}
  \item{model,fixed_method,random_method,Qp_alpha,repeatNum }{value for parmeter of the funciton.}
}
\references{
William S Noble(Nat Biotechnol.2009): How does mutiple testing correction work?

Hedges,L.V. & Vevea,J.L.(1998): Fixed- and random- effects models in meta-analysis.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
'Split':the field separtor of each line for parameter 'genotypeData' and 'affectionData' must be same.
'Split' and 'sep' cannot be same.
}
\seealso{
	\code{\link{permuteGenotype}},
	\code{\link{permuteGenotypeCount}},
	\code{\link{genotypeStat}},
	\code{\link{OR.TradPerm}},
	\code{\link{Armitage.TradPerm}},
   \code{\link{chisq.TradPerm}},
	\code{\link{fisher.TradPerm}},
	\code{\link{meta}},
	\code{\link{meta.MCPerm}},
	\code{\link{VS.Genotype.Hist}},
   \code{\link{VS.Allele.Hist}},
	\code{\link{VS.Hist}},
	\code{\link{PermMeta.LnOR.Hist}},
	\code{\link{PermMeta.LnOR.CDC}},
	\code{\link{PermMeta.LnOR.boxplot}},
	\code{\link{PermMeta.boxplot}},
	\code{\link{PermMeta.Hist}},
	\code{\link{pearson_scatter}},
	\code{\link{Q.TradPerm}},
	\code{\link{I2.TradPerm}}
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
# result=meta.TradPerm(gen,aff,split=" ",sep="/",naString="-",
    # model="allele",method="MH",repeatNum=1000) 
# result
}
\keyword{ meta }
\keyword{ TradPerm }
\keyword{ meta.TradPerm }
