\name{VS.Allele.CDC}
\alias{VS.Allele.CDC}
\title{
separately plot cumulative distribution curve for the return value(allele count) of 'meta.TradPerm' 
and 'meta.MCPerm' for certain study
}
\description{
separately plot cumulative distribution curve for the return value(allele count) of 'meta.TradPerm' 
and 'meta.MCPerm' for certain study.
}
\usage{
VS.Allele.CDC(Trad_case_1, Trad_case_2, Trad_control_1, Trad_control_2, 
    MC_case_1, MC_case_2, MC_control_1, MC_control_2, 
	 Trad_col = "black", MC_col = "red", 
    main = "cumulative distribution curve", 
	 title = c("case_A", "case_a", "control_A", "control_a"), 
	 xlab = "count", ylab = "cumulative probability")
}
\arguments{
  \item{Trad_case_1}{
a numeric vector, simulative allele 1 count for case samples got by TradPerm method for certain study.
}
  \item{Trad_case_2}{
a numeric vector, simulative allele 2 count for case samples got by TradPerm method for certain study.
}
  \item{Trad_control_1}{
a numeric vector, simulative allele 1 count for control samples got by TradPerm method for certain study.
}
  \item{Trad_control_2}{
a numeric vector, simulative allele 2 count for control samples got by TradPerm method for certain study.
}
  \item{MC_case_1}{
a numeric vector, simulative allele 1 count for case samples got by MCPerm method for certain study.
}
  \item{MC_case_2}{
a numeric vector, simulative allele 2 count for case samples got by MCPerm method for certain study.
}
  \item{MC_control_1}{
a numeric vector, simulative allele 1 count for control samples got by MCPerm method for certain study.
}
  \item{MC_control_2}{
a numeric vector, simulative allele 2 count for control samples got by MCPerm method for certain study.
}
  \item{Trad_col}{
the color of cumulative distribution cure for Trad_case_1/Trad_case_2/
Trad_control_1/Trad_control_2. Default value is 'black'.
}
  \item{MC_col}{
the color of cumulative distribution cure for MC_case_1/MC_case_2/
MC_control_1/MC_control_2. Default value is 'red'.
}
  \item{main}{
the main title(on top), default value is "cumulative distribution curve".
}
  \item{title}{
the sub main title for each plot(on top). Default value is a vector with elements: 'case_A', 'case_a', 
'control_A' and 'control_a'.
}
  \item{xlab,ylab}{
X axis label, default value is 'count'.
Y axis label, default value is 'cumulative probability'.
}
}
\details{
Separately plotting cumulative distribution curve for the return value(allele count) of 'meta.TradPerm' 
and 'meta.MCPerm' for certain study is to compare the simulative allele count distribution got by TradPerm 
and MCPerm method whether are same.

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
  \code{\link{VS.QQ}},
  \code{\link{VS.KS}},
  \code{\link{VS.Allele.Hist}},
  \code{\link{VS.Allele.QQ}},
  \code{\link{VS.Genotype.CDC}},
  \code{\link{PermMeta.LnOR.CDC}}
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
# Trad_case_2=2*result1$perm_case_22[12,]+result1$perm_case_12[12,]
# Trad_control_1=2*result1$perm_control_11[12,]+result1$perm_control_12[12,]
# Trad_control_2=2*result1$perm_control_22[12,]+result1$perm_control_12[12,]

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
# MC_case_2=2*result2$perm_case_22[12,]+result2$perm_case_12[12,]
# MC_control_1=2*result2$perm_control_11[12,]+result2$perm_control_12[12,]
# MC_control_2=2*result2$perm_control_22[12,]+result2$perm_control_12[12,]

# VS.Allele.CDC(Trad_case_1,Trad_case_2,Trad_control_1,Trad_control_2,
    # MC_case_1,MC_case_2,MC_control_1,MC_control_2,
    # main="cumulative distribution curve for allele model",
	 # title=c("case_A","case_a","control_A","control_a"))
}
\keyword{ meta.MCPerm }
\keyword{ meta.TradPerm }
\keyword{ Allele }
