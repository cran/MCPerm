\name{PermMeta.LnOR.boxplot}
\alias{PermMeta.LnOR.boxplot}
\title{
boxplot for the return value 'perm_LnOR' or 'perm_VARLnOR' of 'meta.MCPerm' or 'meta.TradPerm'
}
\description{
boxplot for the return value 'perm_LnOR' or 'perm_VARLnOR' of 'meta.MCPerm' or 'meta.TradPerm'.
}
\usage{
PermMeta.LnOR.boxplot(PermMeta, plot = "LnOR", plot_study = "all", order = "no", 
    main = "LnOR,no order", true_value_pch = 3, pch_col = "red", 
	 pos = 3, text_col = "blue", digits = 2)
}
\arguments{
  \item{PermMeta}{
the result of function 'meta.TradPerm' or 'meta.MCPerm'.
}
  \item{plot}{
a character string indicating which return value of function 'meta.TradPerm' or 'meta.MCPerm' to be plot.
The value can be "LnOR"(default), "LnOR_VAR". "LnOR" and "LnOR_VAR" separetly plots the return value 
'perm_LnOR'  and 'perm_VARLnOR'. 
}
  \item{plot_study}{
a numeric vector indicates which study(ies) in meta analysis to be plotted. Default value
is 'all', which indicates all studies in meta analysis to be plotted.
}
  \item{order}{
a character string specifying the boxplot ascending order. The order can be 'LnOR',
'LnOR_VAR', 'VAR_LnOR', 'sample' and 'no'(default). See details.
}
  \item{main}{
the main title (on top), default value is 'LnOR,no order'.
}
  \item{true_value_pch}{ 
the pch(default 3) to mark the observed value in the boxplot of simulation data. 
}
  \item{pch_col}{
the color of pch, default value is 'red'.
}
  \item{pos}{
a position specifier for text of the observed value around the pch.Values of 1, 2, 3(default) and 4, 
respectively indicate positions below, to the left of, above and to the right of the specified 
coordinates.
}
  \item{text_col}{
the color(default 'blue') for the text of the observed value. 
}
  \item{digits}{
integer(default 2) indicating the number of decimal places.
}
}
\details{
Boxplot for the return value('perm_LnOR','perm_VARLnOR') of 'meta.MCPerm' or 'meta.TradPerm'.
And through 'pch' and text to mark the observed value. The boxplot can order by 'LnOR',
'LnOR_VAR', 'VAR_LnOR', 'sample' and 'no'. 

'LnOR' indicates the boxplot order by the log odd ratio calculated by each study 
observed genotype count.  'LnOR_VAR' indicates the boxplot order by the variance of 
log odd ratio calculated by formule 1/ai+1/bi+1/ci+1/di(ai/bi/ci/di separetly are 
case/control two alleles count).  'VAR_LnOR' indicates the variance of simulation 
log odd ratios.  'sample' indicates the boxplot order by the sample size of each study.
'no' indicates that the boxplot order is same as the order of parameter 'plot_study'.

MCPerm details see \code{\link{chisq.MCPerm}}.
TradPerm details see \code{\link{chisq.TradPerm}}.
}
\value{
  \item{plot }{character value of parameter 'plot'.}
  \item{plot_num }{the number of the selected studies to be plotted.}
  \item{plot_order }{character value of parmeter 'order'.}
  \item{order_index }{a numeric vector of the order of studies.}
  \item{true_LnOR }{a numeric vector of log odd ratio calculated by observed value of each 
   study. And the order is corresponding to the 'order_index'.}
  \item{LnOR_VAR }{a numeric vector of the log odd ratio variance(1/ai+1/bi+1/ci+1/di) 
     calculated by observed value of each study. And the order is corresponding to the 'order_index'.  }
  \item{VAR_LnOR }{a numeric vector of the variance of the simulated log odd ratios got by permutation.
     And the order is corresponding to the 'order_index'.}
  \item{sample }{a numeric vector of the sample size of the plotted studies, and the order is 
    corresponding to the 'order_index'.}
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
## set working directory to save the plots.
# setwd("D:\")
# pdf("PermMeta.LnOR.boxplot.pdf",height=6,width=9)
# par(mfrow=c(2,3))
# PermMeta.LnOR.boxplot(result,plot="LnOR",plot_study="all",order="no")
# PermMeta.LnOR.boxplot(result,plot="LnOR",plot_study="all",order="LnOR",
    # main="LnOR, order by LnOR")
# PermMeta.LnOR.boxplot(result,plot="LnOR",plot_study="all",order="LnOR_VAR",
    # main="LnOR, order by LnOR_VAR")
# PermMeta.LnOR.boxplot(result,plot="LnOR",plot_study="all",order="VAR_LnOR",
    # main="LnOR, order by VAR_LnOR")
# PermMeta.LnOR.boxplot(result,plot="LnOR",plot_study="all",order="sample",
    # main="LnOR, order by sample")

# PermMeta.LnOR.boxplot(result,plot="LnOR_VAR",plot_study=c(1,4,14,12),order="no",
    # main="LnOR_VAR, no order")
# PermMeta.LnOR.boxplot(result,plot="LnOR_VAR",plot_study=c(1,4,14,12),order="LnOR",
    # main="LnOR_VAR, order by LnOR")
# PermMeta.LnOR.boxplot(result,plot="LnOR_VAR",plot_study=c(1,4,14,12),order="LnOR_VAR",
    # main="LnOR_VAR, order by LnOR_VAR")
# PermMeta.LnOR.boxplot(result,plot="LnOR_VAR",plot_study=c(1,4,14,12),order="VAR_LnOR",
    # main="LnOR_VAR, order by VAR_LnOR")
# PermMeta.LnOR.boxplot(result,plot="LnOR_VAR",plot_study=c(1,4,14,12),order="sample",
    # main="LnOR_VAR, order by sample")
# dev.off()
}

\keyword{ PermMeta }
\keyword{ meta.MCPerm }
\keyword{ meta.TradPerm }