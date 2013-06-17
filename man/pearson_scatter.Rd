\name{pearson_scatter}
\alias{pearson_scatter}
\title{
scatter plot and calculate Pearson correlation coefficient for paired data
}
\description{
scatter plot and calculate Pearson correlation coefficient for paired data.
}
\usage{
pearson_scatter(Trad_data, MC_data, scatter_col = "gray28", line_col = "black", 
    title = NULL, xlab = "TradPerm P-value", ylab = "MCPerm P-value")
}
\arguments{
  \item{Trad_data}{
numeric vector, e.g. the result('perm_Qp'/'perm_I2'/'perm_p') of function 'meta.TradPerm'.
}
  \item{MC_data}{
numeric vector, e.g. the result('perm_Qp'/'perm_I2'/'perm_p') of function 'meta.MCPerm'.
}
  \item{scatter_col}{
the color(default 'gray28') of the scatter points.
}
  \item{line_col}{
the color(default 'black') of the line x=y.
}
  \item{title}{
The main title (on top).
}
  \item{xlab,ylab}{
X axis label, default value is 'TradPerm P-value'.
Y axis label, default value is 'MCPerm P-value'.
}
}
\details{
Scatter plot and Pearson correlation coefficient(two.sided) for 'perm_Qp'/'perm_I2'/'perm_p' of
'meta.TradPerm' and 'meta.MCPerm' are to test the consistency between them.
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
  \code{\link{PermMeta.Hist}},
  \code{\link{PermMeta.boxplot}}
}
\examples{
# Trad=read.table("Trad_result.txt",sep=" ",header=FALSE)
# MC=read.table("MC_result.txt",sep=" ",header=FALSE)
# par(mfrow=c(3,1))
# pearson_scatter(as.numeric(Trad[,4]),as.numeric(MC[,4]),
    # title="Q p_value ",
    # xlab="TradPerm Qp_value",
    # ylab="MCPerm Qp_value")
# pearson_scatter(as.numeric(Trad[,6]),as.numeric(MC[,6]),
    # title="I2 p_value",
    # xlab="TradPerm I2p_value",
    # ylab="MCPerm I2p_value")	 
# pearson_scatter(as.numeric(Trad[,9]),as.numeric(MC[,9]),
    # title="p_value",
    # xlab="TradPerm p_value",
    # ylab="MCPerm p_value")
}
\keyword{ meta.MCPerm }
\keyword{ meta.TradPerm }
