\name{permuteGenotype}
\alias{permuteGenotype}
\title{
Permute the elements of genotype data
}
\description{
Permute the elements of genotype data.
}
\usage{
permuteGenotype(dataLine, fromCol)
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
Permute the elements of genotype data.
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
   \code{\link{Armitage.TradPerm}},
	\code{\link{OR.TradPerm}},
	\code{\link{chisq.TradPerm}},
	\code{\link{fisher.TradPerm}},
	\code{\link{meta.TradPerm}},
	\code{\link{permuteGenotypeCount}},
	\code{\link{chisq.MCPerm}}
}
\examples{
# dataLine1=matrix(c("rs12","AA","AG","GG","CG","AA"),nrow=1)
# permuteGenotype(dataLine1,fromCol=2)
# dataLine2=matrix(c("rs12","AA","AG","GG","CG","AA",
   # "rs12","AA","AG","GG","CG","AA"),nrow=2,byrow=TRUE)
# permuteGenotype(dataLine2,fromCol=2)
}
\keyword{ permutation }
\keyword{ genotype }
