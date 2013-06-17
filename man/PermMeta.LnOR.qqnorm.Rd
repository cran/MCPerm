\name{PermMeta.LnOR.qqnorm}
\alias{PermMeta.LnOR.qqnorm}
\title{
qqnorm plot for the return value 'perm_LnOR' of 'meta.MCPerm' or 'meta.TradPerm'
}
\description{
qqnorm plot for the return value 'perm_LnOR' of 'meta.MCPerm' or 'meta.TradPerm'.
}
\usage{
PermMeta.LnOR.qqnorm(PermMeta, plot_study = "all", nrow = 2, ncol = 2, 
    main = "qqnorm plot for LnOR", title = NULL, 
	 xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", 
	 scatter_col = "black", line_col = "red")
}
\arguments{
  \item{PermMeta}{
the result of function 'meta.TradPerm' or 'meta.MCPerm'.
}
  \item{plot_study}{
a numeric vector indicates which study(ies) in meta analysis to be plotted. Default value is 
'all', which indicates all studies in meta analysis to be plotted.
}
  \item{nrow,ncol}{
positive integer, divides the device up into 'nrow' rows and 'ncol' columns. Default value is 2.
}
  \item{main}{
the main title (on top), default value is "qqnorm plot for LnOR".
}
  \item{title}{
the sub main title for each plotted study(on top).
}
  \item{xlab,ylab}{
X axis label, default value is "Theoretical Quantiles".
Y axis label, default value is "Sample Quantiles".
}
  \item{scatter_col}{
the color of the scatter points. Default value is 'black'.
}
  \item{line_col}{
the color of the line which passes through the normal distribution probs quantiles, the first and 
third quartiles. Default value is 'red'.
}
}
\details{
Plotting a normal QQ plot for simulative log odd ratios is to test that simulative data is whether 
fit normal distribution, Snd plot a line which passes through the normal distribution probs quantiles,
the first and third quartiles.

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
  \code{\link{VS.QQ}},
  \code{\link{VS.KS}},
  \code{\link{VS.Genotype.QQ}},
  \code{\link{VS.Allele.QQ}},
  \code{\link{PermMeta.LnOR.Hist}},
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
# PermMeta.LnOR.qqnorm(result,plot_study=c(1,4,9,15),nrow=2,ncol=2)
}
\keyword{ PermMeta }
