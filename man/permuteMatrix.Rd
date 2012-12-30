\name{permuteMatrix}
\alias{permuteMatrix}
\title{
 Permutate Matrix with two rows and row/columns sums fixed
}
\description{
 Fixed the rows and columns sums,then fill matrix.
}
\usage{
permuteMatrix(rawMatrix)
}
\arguments{
  \item{rawMatrix}{
a matrix with two rows.
}
}
\details{
For a matrix with two rows, fix it rows and columns sums,then fill matrix by generating random number  
for the hypergeometric distribution. 
}
\value{
The return value is a matrix with the same dimension with the inputting matrix.
}
\references{
Edgington. E.S. (1995) Randomization tests, 3rd ed. New York: Marcel-Dekker
}
\author{
Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
The inputting matrix must be two rows.
}
\seealso{
  \code{\link{genotypeStat}},
  \code{\link{permuteData}},   
  \code{\link{OR}},   
  \code{\link{Armitage}}, 
  \code{\link{tradPerm}},
  \code{\link{MCPerm}},
  \code{\link{rhyper}}
}
\examples{
rawMatrix=matrix(c(24,7,8,34),nrow=2)
newMatrix=permuteMatrix(rawMatrix)
}
\keyword{permute }
\keyword{Matrix}
