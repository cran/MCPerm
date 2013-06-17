\name{PermMeta.LnOR.CDC}
\alias{PermMeta.LnOR.CDC}
\title{
cumulative distribution curve for the return value 'perm_LnOR' of 'meta.MCPerm' or
'meta.TradPerm'
}
\description{
cumulative distribution curve for the return value 'perm_LnOR' of 'meta.MCPerm' or
'meta.TradPerm'.
}
\usage{
PermMeta.LnOR.CDC(PermMeta, plot_study = "all", nrow = 2, ncol = 3, 
    PermMeta.LnOR_pch = 4, PermMeta.LnOR_col = "black", 
	 LnOR_VAR_pch = 18, LnOR_VAR_col = "blue", 
	 VAR_LnOR_pch = 18, VAR_LnOR_col = "red", 
	 main = "cumulative distribution curve for LnOR", title = NULL, 
	 xlab = "LnOR", ylab = "cumulative probability", digits = 3)
}
\arguments{
  \item{PermMeta}{
the result of function 'meta.TradPerm' or 'meta.MCPerm'.
}
  \item{plot_study}{
a numeric vector indicates which study(ies) in meta analysis to be plotted. Default value
is 'all', which indicates all studies in meta analysis to be plotted.
}
  \item{nrow,ncol }{
positive integer, divides the device up into 'nrow'(default is 2) rows and 'ncol'(default is 3) columns. 
}
  \item{PermMeta.LnOR_pch,PermMeta.LnOR_col }{
the pch(default 4) and the color(default 'black') of pch are for the cumulative distribution curve of the 
return value 'perm_LnOR' of certain study. 
}
  \item{LnOR_VAR_pch,LnOR_VAR_col }{
the pch(default 18) and the color*(default 'blue') of pch are for the cumulative distribution curve of the normal distrition
with mean=0 and variance=1/ai+1/bi+1/ci+1/di.
}
  \item{VAR_LnOR_pch,VAR_LnOR_col }{
the pch(default 18) and the color(default 'red') of pch are for the cumulative distribution curve of the normal distrition
with mean=0 and variance is the variance of simulation log odd ratios.
}
  \item{main}{
the main title(on top), default value is "cumulative distribution curve for LnOR". 
}
  \item{title}{
the sub main title for each plotted study(on top).
}
  \item{xlab,ylab}{
X axis label, default value is 'LnOR'.
Y axis label, default value is 'cumulative probability'.
}
  \item{digits}{
integer(default 3) indicating the number of decimal places.
}
}
\details{
Plot three cumulative distribution cures(abbreviation:CDC): 1) CDC for simulative log odd ratios;
2) CDC for normal distribution with mean=0 and var=variance of observed log odd ratio of
certain study(1/ai+1/bi+1/ci+1/di); 3) CDC for normal distribution with mean=0 and var=variance 
of simulative log odd ratios. The symbol---'perm_lnOR','pnorm_LnOR_VAR','pnorm_VAR_LnOR' in 
the topright legend separately indicated the first, second, third cure. Through three CDC 
compared, observe than the thrid cure is more corresponding to the first cure when sample 
size is smaller. 

The symbol in the bottomright legend: 'LnOR' indicates the log odd ratio of observed data
for the study; 'sample' indicates the sample size of the study; 'LnOR_VAR' indicates the 
variance of second cure; 'VAR_LnOR' indicates the variance of third cure. 

MCPerm details see \code{\link{chisq.MCPerm}}.
TradPerm details see \code{\link{chisq.TradPerm}}.
}
\value{
  \item{plot_study }{the value of paramter 'plot_study'.}
  \item{LnOR }{the numeric vector of log odd ratio of the observed data for the plotted studies. }
  \item{sample }{the numeric vector of sample size of the plotted studies.}
  \item{LnOR_VAR }{the numeric vector of variance of the sencond cure for the plotted studies. }
  \item{VAR_LnOR }{the numeric vector of variance of the third cure for the plotted studies.}
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
  \code{\link{meta.MCPerm}},
  \code{\link{meta.TradPerm}},
  \code{\link{chisq.MCPerm}},
  \code{\link{chisq.TradPerm}},
  \code{\link{VS.CDC}},
  \code{\link{VS.KS}},
  \code{\link{VS.Genotype.CDC}},
  \code{\link{VS.Allele.CDC}},
  \code{\link{PermMeta.LnOR.Hist}},
  \code{\link{PermMeta.LnOR.qqnorm}},
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
# PermMeta.LnOR.CDC(result,plot_study=c(3,5,21,7,12,9),nrow=2,ncol=3)
}
\keyword{ PermMeta }
