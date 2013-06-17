\name{PermMeta.Hist}
\alias{PermMeta.Hist}
\title{
histplot for the result of 'meta.MCPerm' or 'meta.TradPerm'
}
\description{
histplot for the result of 'meta.MCPerm' or 'meta.TradPerm'.
}
\usage{
PermMeta.Hist(PermMeta, plot = "Qp", fill_col = NULL, border_col = "black", 
    arrows_col = "red", main = "Hist plot for heterogeneity Q p_vlaue", 
	 xlab = "Q p_value", ylab = "Density", digits = 3)
}
\arguments{
  \item{PermMeta}{
the result of function 'meta.TradPerm' or 'meta.MCPerm'.
}
  \item{plot}{
a character string indicating which return value of function 'meta.TradPerm' or 'meta.MCPerm' to be plot.
The value can be "Qp"(default), "I2", "merged_LnOR", "merged_LnOR_VAR" or "merged_LnOR_p". And the value must be simulation data.
'Qp', "I2", "merged_LnOR", "merged_LnOR_VAR" and "merged_LnOR_p" separetly plots the return value 'perm_Qp', 'perm_I2',
'perm_merged_LnOR', 'perm_merged_VARLnOR', 'perm_merged_p'.
}
  \item{fill_col}{
the filled color(default NULL) for the body of histplot.
}
  \item{border_col}{
the color(default 'black') for the border of histplot.
}
  \item{arrows_col}{
the color(default 'red') of arrows which mark the place of the observed value.
}
  \item{main}{
the main title (on top), default value is "Hist plot for heterogeneity Q p_vlaue".
}
  \item{xlab,ylab}{
X axis label, default value is 'Q p_value'.
Y axis label, default value is 'Density'.
}
  \item{digits}{
integer(default 3) indicating the number of decimal places.
}
}
\details{
Histplot for the return value('perm_Qp','perm_I2','perm_merged_LnOR','perm_merged_VARLnOR',
'perm_merged_p') of 'meta.MCPerm' or 'meta.TradPerm'. And through arrows and legend to mark the 
observed value.

The symbols in the legend: 'Q_stat' is the Q statistic for meta data heterogeneity; 
'Q_p' is the p value of Q value(chi square distribution,the number of studies in 
meta analysis minus one is degree of freedom of Q value.); 'p.corrected' is the corrected
p value by permutation; 'I2_stat' is the statistic I2(calculated by formula max(Q-d.f./Q, 0))
for meta data heterogeneity; 'merged_LnOR' is the merged log odd ratio of observed data; 
'merged_LnOR_VAR' is the variance of log odd ratio for observed data; 'merged_LnOR_p' is 
the p value of log odd ratio of observed data which obey normal distribution. 

MCPerm details see \code{\link{chisq.MCPerm}}.
TradPerm details see \code{\link{chisq.TradPerm}}.
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
  \code{\link{PermMeta.LnOR.Hist}},
  \code{\link{PermMeta.LnOR.CDC}},
  \code{\link{PermMeta.LnOR.boxplot}},
  \code{\link{PermMeta.LnOR.qqnorm}},
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
## set working directory to save the plots.
# setwd("D:\")
# pdf("PermMeta.Hist.pdf",height=6,width=6)
# PermMeta.Hist(result,plot="Qp",fill_col=NULL,border_col='black',
    # arrows_col='red',main="Hist plot for heterogeneity Q p_vlaue",xlab="Q p_value")
# PermMeta.Hist(result,plot="I2",fill_col=NULL,border_col='black',
    # arrows_col='red',main="Hist plot for heterogeneity I2",xlab="I2")
# PermMeta.Hist(result,plot="merged_LnOR",fill_col=NULL,border_col='black',
    # arrows_col='red',main="Hist plot for merged_LnOR",xlab="merged_LnOR")
# PermMeta.Hist(result,plot="merged_LnOR_VAR",fill_col=NULL,border_col='black',
    # arrows_col='red',main="Hist plot for merged_LnOR_VAR",xlab="merged_LnOR_Variance")
# PermMeta.Hist(result,plot="merged_LnOR_p",fill_col=NULL,border_col='black',
    # arrows_col='red',main="Hist plot for merged_LnOR_p",xlab="merged_LnOR_p.value")
# dev.off()
}
\keyword{ PermMeta }
\keyword{ MCPerm }