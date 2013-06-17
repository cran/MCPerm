\name{meta}
\alias{meta}
\title{
Meta analysis
}
\description{
Meta analysis.
}
\usage{
meta(model, fixed_method, random_method, Qp_alpha, case_11, case_12, case_22, 
    control_11, control_12, control_22, label = NULL, dataset = NULL)
}
\arguments{
  \item{model}{
a character string indicating the type of model("allele","dominant" or "recessive") supplied to the data.
The risk allele(see details) is marked as allele1. The allele model indicates allele1 versus allele2, the 
dominant model indicates <allele1/allele1> + <allele1/allele2> versus <allele2/allele2>,
the recessive model indicates <allele1/allele1> versus <allele1/allele2> + <allele2/allele2>.
}
  \item{fixed_method}{
a character string indicating the method('Inverse','MH' or 'Peto') to fit fixed-effects model.
}
  \item{random_method}{
a character string indicating the method to fit random-effects model: "HE","DL",
"SJ","ML","REML" or "EB".
}
  \item{Qp_alpha}{
the threshold value(0-1) to refuse the null hypothesis that all studies were evaluating the same effect.
Not rejecting the above hypothesis usually leads a meta-analysis to adopt a fixed-effects model.
The fixed-effects model assumes that the estimated effect sizes only differ by the sampling error. 
In contrast, if a significant Q-statistic P<Qp_alpha indicates heterogeneity across studies, a 
random-effects model should be adopted.
}
  \item{case_11}{
non-negative integer vector, the frequency of genotype <allele1/allele1> in case samples.
}
  \item{case_12}{
non-negative integer vector, the frequency of genotype <allele1/allele2> in case samples.
}
  \item{case_22}{
non-negative integer vector, the frequency of genotype <allele2/allele2> in case samples.
}
  \item{control_11}{
non-negative integer vector, the frequency of genotype <allele1/allele1> in control samples.
}
  \item{control_12}{
non-negative integer vector, the frequency of genotype <allele1/allele2> in control samples.
}
  \item{control_22}{
non-negative integer vector, the frequency of genotype <allele2/allele2> in control samples.
}
   \item{label}{
    character string vector, labels for the studies to show in the forest plot. When label=NULL(default),
	the label for the ith study will be "study i". Also see parameter 'dataset'.
}
  \item{dataset}{
data frame containing the data of genotype frequency(e.g. 'case_11','case_12') and the information 
of the studies(e.g. author,year). When dataset=NULL(default), the value of parameters 'case_11', 'case_12', 
'case_22', 'control_11', 'control_12', 'control_22', 'label' must be vector with same length(>1).
When dataset is not NULL, the value of parameters 'case_11', 'case_12', 'case_22', 'control_11', 'control_12',
'control_22' must be a non-negative integer, that specify the column of parameter needed value in 'dataset';
and lable can be a non-negative integer or vector that will be connected by " " as the labels for studies.
}
}
\details{
Allele 1 and allele 2 to each study have OR values. The risk allele is the allele which 
the number of studies which OR>1 more than half of the number of all studies.

Meta-analysis refers to methods focused on contrasting and combining results from different studies, 
in the hope of identifying patterns among study results, sources of disagreement among those results, 
or other interesting relationships that may come to light in the context of multiple studies. In its 
simplest form, this is normally by identification of a common measure of effect size, of which a weighted 
average might be the output of a meta-analysis. The weighting might be related to sample sizes within 
the individual studies. More generally there are other differences between the studies that need to be 
allowed for, but the general aim of a meta-analysis is to more powerfully estimate the true effect 
size as opposed to a less precise effect size derived in a single study under a given single set 
of assumptions and conditions.
}
\value{
The return value detail see package metafor. See package:metafor.
If the return result is fixed-effect model, see metafor:rma(method='FE'), metafor:rma.mh, metafor:rma.peto.
If the return result is random-effect model, see metafor:rma.
   
   \item{b }{combined log odds ratio.}
   \item{se }{standard errors of the combined log odds ratio.}
   \item{zval }{test statistics of the combined log odds ratio.}
   \item{pval }{p-values for the test statistics.}
   \item{ci.lb }{lower bound of the confidence intervals for the combined log odds ratio.}
   \item{ci.ub }{upper bound of the confidence intervals for the combined log odds ratio.}
   \item{QE }{test statistic for the test of heterogeneity.}
   \item{QEp }{p-value for the test of heterogeneity.}
}
\references{
William S Noble(Nat Biotechnol.2009): How does mutiple testing correction work?

Hedges,L.V. & Vevea,J.L.(1998): Fixed- and random- effects models in meta-analysis.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
	\code{\link{meta.TradPerm}},
	\code{\link{meta.MCPerm}},
	\code{\link{permuteGenotype}},
	\code{\link{permuteGenotypeCount}},
	\code{\link{genotypeStat}},
	\code{\link{OR.TradPerm}},
	\code{\link{Armitage.TradPerm}},
   \code{\link{chisq.TradPerm}},
	\code{\link{fisher.TradPerm}},
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
# case_11=c(1,2,3,4,7)
# case_12=c(1,2,2,2,4)
# case_22=c(1,2,3,4,5)
# control_11=c(1,4,5,6,3)
# control_12=c(3,4,5,6,6)
# control_22=c(3,6,7,8,5)
# result1=meta("allele","Inverse","DL",0.05,case_11,case_12,case_22,
    # control_11,control_12,control_22)
# result2=meta("dominant","MH","DL",0.05,case_11,case_12,case_22,
    # control_11,control_12,control_22)
# result3=meta("recessive","Peto","DL",0.05,case_11,case_12,case_22,
    # control_11,control_12,control_22)
}
\keyword{ meta }
