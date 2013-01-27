\name{permuteData}
\alias{permuteData}
\title{
Permute the elements of a matrix
}
\description{
Permute the elements of a matrix.
}
\usage{
permuteData(dataLine, fromCol)
}
\arguments{
  \item{dataLine}{
A matrix with one row.
}
  \item{fromCol}{
A positive integer, the start column to permute.
}
}
\details{
Permute the elements of a matrix.
}
\value{
The return value is a matrix with one row and the elements has been permutated.
}
\references{
William S Noble(Nat Biotechnol.2009): How does multiple testing correction work? 

Edgington. E.S.(1995): Randomization tests, 3rd ed.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
    \code{\link{Armitage}},
	\code{\link{Armitage.TradPerm}},
    \code{\link{Armitage.MCPerm}},
    \code{\link{OR}},
	\code{\link{OR.TradPerm}},
	\code{\link{OR.MCPerm}},
	\code{\link{permuteGenotype}},
	\code{\link{genotypeStat}},
    \code{\link{chisq.TradPerm}},
	\code{\link{chisq.MCPerm}},
	\code{\link{fisher.TradPerm}},
	\code{\link{fisher.MCPerm}},
	\code{\link{rhyper}},
	\code{\link{chisq.test}},
	\code{\link{fisher.test}}
}
\examples{
dataLine1=matrix(c("rs12","AA","AG","GG","CG","AA"),nrow=1)
permuteData(dataLine1,fromCol=2)
dataLine2=matrix(c("rs12","AA","AG","GG","CG","AA","rs12","AA","AG","GG","CG","AA"),nrow=2,byrow=TRUE)
permuteData(dataLine2,fromCol=2)
}
\keyword{ permutation }
