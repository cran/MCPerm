\name{PermMeta.LnOR.Hist}
\alias{PermMeta.LnOR.Hist}
\title{
histplot for the return value 'perm_LnOR' or 'perm_VARLnOR' of 'meta.MCPerm' or
'meta.TradPerm'
}
\description{
histplot for the return value 'perm_LnOR' or 'perm_VARLnOR' of 'meta.MCPerm' or
'meta.TradPerm'
}
\usage{
PermMeta.LnOR.Hist(PermMeta, plot = "LnOR", plot_study = "all", nrow = 2, ncol = 2, 
    main = "Background distribution for LnOR", title = NULL, 
	 xlab = "LnOR", hist_border_col = "black", arrows_col = "red", digits = 3)
}
\arguments{
  \item{PermMeta}{
the result of function 'meta.TradPerm' or 'meta.MCPerm'.
}
  \item{plot}{
a character string indicating which return value of function 'meta.TradPerm' or 
'meta.MCPerm' to be plot. The value can be "LnOR", "LnOR_VAR". "LnOR" and 
"LnOR_VAR" separetly plots the return value 'perm_LnOR'  and 'perm_VARLnOR'.
Default value is 'LnOR'.
}
  \item{plot_study}{
a numeric vector indicates which study(ies) in meta analysis to be plotted. Default value
is 'all', which indicates all studies in meta analysis to be plotted.
}
  \item{nrow,ncol}{
positive integer, divides the device up into 'nrow' rows and 'ncol' columns. 
Default value is 2.
}
  \item{main}{
the main title (on top), default value is "Background distribution for LnOR".
}
  \item{title}{
the sub main title for each plotted study(on top).
}
  \item{xlab}{
X axis label, default value is 'LnOR'.
}
  \item{hist_border_col}{
the color for the border of histplot. Default value is 'black'.
}
  \item{arrows_col}{
the col of arrows which mark the place of the observed value. Default value is 'red'.
}
  \item{digits}{
integer indicating the number of decimal places. Default value is 3.
}
}
\details{
Histplot for the return value('perm_LnOR','perm_VARLnOR') of 'meta.MCPerm' or 'meta.TradPerm'. 
And through arrows and legend to mark the observed value.

The symbol in the legend: 'LnOR' indicates the log odd ratio of observed value for the study.  
'LnOR_VAR' indicates the variance of log odd ratio calculated by formule 1/ai+1/bi+1/ci+1/di
(ai/bi/ci/di separetly are case/control two alleles count). 'sample' indicates the sample size
of the study. 'mean' indicates the mean value of simulative log odd ratios. 'var' indicates the
variance of simulative log odd ratios. 'p' is the percentage of more than the observed value.

MCPerm details see \code{\link{chisq.MCPerm}}.
TradPerm details see \code{\link{chisq.TradPerm}}.
}
\value{
  \item{plot_study }{the value of paramter 'plot_study'.}
  \item{LnOR }{the numeric vector of log odd ratio of the observed data for the plotted studies. }
  \item{sample }{the numeric vector of sample size of the plotted studies.}
  \item{LnOR_VAR }{the numeric vector of variance of the log odd ratio of obaserved value for the
    plotted studies. }
  \item{VAR_LnOR }{the numeric vector of variance of simulative log odd ratios for the plotted studies.}
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
  \code{\link{meta.MCPerm}},
  \code{\link{meta.TradPerm}},
  \code{\link{chisq.MCPerm}},
  \code{\link{chisq.TradPerm}},
  \code{\link{VS.Hist}},
  \code{\link{VS.KS}},
  \code{\link{VS.Genotype.Hist}},
  \code{\link{VS.Allele.Hist}},
  \code{\link{PermMeta.LnOR.qqnorm}},
  \code{\link{PermMeta.LnOR.CDC}},
  \code{\link{PermMeta.LnOR.boxplot}},
  \code{\link{PermMeta.Hist}},
  \code{\link{PermMeta.boxplot}}
}
\examples{
## import data
# data(MetaGenotypeCount)
## delete first line
# temp=MetaGenotypeCount[-1,];
# result=meta.MCPerm(case_11=as.numeric(temp[,14]),case_12=as.numeric(temp[,16]),
	 # case_22=as.numeric(temp[,18]),control_11=as.numeric(temp[,15]),
	 # control_12=as.numeric(temp[,17]),control_22=as.numeric(temp[,19]),
	 # model="allele",fixed_method="MH",random_method="DL",repeatNum=1000)
# PermMeta.LnOR.Hist(result,plot="LnOR",plot_study="all",nrow=2,ncol=2,
    # title=NULL,xlab="LnOR")
# PermMeta.LnOR.Hist(result,plot="LnOR_VAR",plot_study=c(3,5,12,21),nrow=2,ncol=2,
    # main="Background distribution for LnOR_VAR",title=NULL,xlab="LnOR_VAR")
}
\keyword{ PermMeta }
