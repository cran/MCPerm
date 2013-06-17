\name{VS.QQ}
\alias{VS.QQ}
\title{
plot quantile-quantile plot for the return value of 'meta.TradPerm' 
and 'meta.MCPerm' for certain study or meta analysis.
}
\description{
plot quantile-quantile plot for the return value of 'meta.TradPerm' 
and 'meta.MCPerm' for certain study or meta analysis.
}
\usage{
VS.QQ(Trad_data, MC_data, scatter_col = "black", line_col = "black", title = "QQ plot", 
    xlab = "Quantile for TradPerm data)", ylab = "Quantile for MCPerm data")
}
\arguments{
  \item{Trad_data}{
the return value of function 'meta.TradPerm', e.g. 'perm_case_11' of certain stuy, 'perm_Qp', 'perm_p' etc.
}
  \item{MC_data}{
the return value of function 'meta.MCPerm', e.g. 'perm_case_11' of certain stuy, 'perm_Qp', 'perm_p' etc.
}
  \item{scatter_col}{
the color for scatter points of quantile-quantile plot, default value is 'black'.
}
  \item{line_col}{
the color of line which passes through the sample distribution probs quantiles, the first and 
third quartiles. Default value is 'black'.
}
  \item{title}{
the main title(on top), default value is 'QQ plot'.
}
  \item{xlab,ylab}{
X axis label, default is "Quantile for TradPerm data)". Y axis label, default is "Quantile for MCPerm data".
}
}
\details{
Plotting quantile-quantile plot for the return value(e.g. 'perm_case_11' of certain stuy,
'perm_Qp', 'perm_p' etc) of 'meta.TradPerm' and 'meta.MCPerm' is to compare the simulative data distribution 
got by TradPerm and MCPerm method whether are same.

MCPerm details see \code{\link{chisq.MCPerm}}.
TradPerm details see \code{\link{chisq.TradPerm}}.
}
\references{
William S Noble(Nat Biotechnol.2009): How does multiple testing correction work?

Edgington. E.S.(1995): Randomization tests, 3rd ed.
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
  \code{\link{VS.CDC}},
  \code{\link{VS.KS}},
  \code{\link{VS.Allele.QQ}},
  \code{\link{VS.Genotype.QQ}}
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
# result1=meta.TradPerm(gen,aff,split=" ",sep="/",naString="-",
    # model="allele",method="MH",repeatNum=1000) 
# result1
## plot study 12
# Trad_case_1=2*result1$perm_case_11[12,]+result1$perm_case_12[12,]

## import data
# data(MetaGenotypeCount)
## delete the first line which is the names for columns.
# temp=MetaGenotypeCount[-1,,drop=FALSE]
# result=meta.MCPerm(case_11=as.numeric(temp[,14]),case_12=as.numeric(temp[,16]),
	 # case_22=as.numeric(temp[,18]),control_11=as.numeric(temp[,15]),
	 # control_12=as.numeric(temp[,17]),control_22=as.numeric(temp[,19]),
	 # model="allele",method="MH",repeatNum=100000)
# result2
## plot study 12
# MC_case_1=2*result2$perm_case_11[12,]+result2$perm_case_12[12,]

# VS.QQ(Trad_case_1,MC_case_1,title="cumulative distribution cure for case_1")
# VS.QQ(result1$perm_Qp,result2$perm_Qp,title="cumulative distribution cure for Qp")
# VS.QQ(result1$perm_p,result2$perm_p,title="cumulative distribution cure for p")

}
\keyword{ meta.TradPerm }
\keyword{ meta.MCPerm }
