\name{VS.KS}
\alias{VS.KS}
\title{
Kolmogorov-Smirnov test for the return value of 'meta.TradPerm' and 'meta.MCPerm'
}
\description{
Kolmogorov-Smirnov test for the return value of 'meta.TradPerm' and 'meta.MCPerm'.
}
\usage{
VS.KS(Trad_data, MC_data, scatter_alpha = 0.01, line_alpha = 0.001, 
    scatter_col = "black", line_col = "red", 
	 xlab = NULL, ylab = "KS test p_value", 
	 title = "KS (Kolmogorov-Smirnov) test")
}
\arguments{
  \item{Trad_data}{
a matrix with more than one line, each line is the return value of function 'meta.TradPerm'(e.g. 'perm_case_11' of 
certain study, 'perm_Qp', 'perm_p' etc).
}
  \item{MC_data}{
a matrix with more than one line, each line is the return value of function 'meta.MCPerm'(e.g. 'perm_case_11' of 
certain study, 'perm_Qp', 'perm_p' etc).
}
  \item{scatter_alpha}{
numeric indicates k-s test p.value less 'scatter_alpha' to be plotted. Default value is 0.01.
}
  \item{line_alpha}{
numeric indicates the threshold value of k-s test significance. Default value is 0.001.
}
  \item{scatter_col}{
the color of sctter points which to be plotted. Default value is 'black'.
}
  \item{line_col}{
the color of the line. Default value is 'red'.
}
  \item{xlab,ylab}{
X axis label. Y axis label, default value is 'KS test p_value'.
}
  \item{title}{
the main title(on top),default value is "KS (Kolmogorov-Smirnov) test".
}
}
\details{
Kolmogorov-Smirnov test for the return value(e.g. 'perm_case_11', 'perm_Qp' etc.) of 'meta.TradPerm' and 
'meta.MCPerm' to test the corresponding form the whole. Just plotting p.value less the 'scatter_alpha' will
more clear to see the number of no-corresponding data.
}
\value{
  \item{KS_p }{a numeric vector, p.value of K-S test.}
}
\references{
The two-sided one-sample distribution comes via Marsaglia, Tsang and Wang (2003). 
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
Parameter 'Trad_data' and 'MC_data' are matrix with more than one line.
}
\seealso{
  \code{\link{meta.MCPerm}},
  \code{\link{meta.TradPerm}},
  \code{\link{chisq.MCPerm}},
  \code{\link{chisq.TradPerm}},
  \code{\link{VS.Hist}},
  \code{\link{VS.CDC}},
  \code{\link{VS.QQ}}
  \code{\link{VS.Allele.Hist}},
  \code{\link{VS.Genotype.Hist}}
}
\examples{
## write the return value 'perm_case_11','perm_Qp' or 'perm_p' of function 'meta.TradPerm' 
  # to file 'Trad_data.txt'
## write the return value 'perm_case_11','perm_Qp' or 'perm_p' of function 'meta.MCPerm' 
  # to file 'MC_data.txt'
## when all meta analysis run over, read the file to K-S test
## if the line of file is separated by '\t'
# Trad=read.table("Trad_data.txt",sep="\t",header=FALSE)
# MC=read.table("MC_data.txt",sep="\t",header=FALSE)

# VS.KS(Trad, MC)
}
\keyword{ meta.MCPerm }
\keyword{ meta.TradPerm }
