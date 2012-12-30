\name{permuteData}
\alias{permuteData}
\title{
Permute Data
}
\description{
Permute Data.
}
\usage{
permuteData(dataLine, fromCol)
}
\arguments{
  \item{dataLine}{
matrix or data.frame with one row need to be permuted.
}
  \item{fromCol}{
integer,the start column to permute.
}
}
\details{
Permute Data. 
}
\value{
The return value is a data.frame with one row.
}
\references{
Edgington. E.S. (1995) Randomization tests, 3rd ed. New York: Marcel-Dekker
}
\author{
Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
  \code{\link{genotypeStat}},
  \code{\link{permuteMatrix}},   
  \code{\link{OR}}, 
  \code{\link{Armitage}},  
  \code{\link{tradPerm}},
  \code{\link{MCPerm}}
}
\examples{
dataLine=matrix(c("rs12","1","2","3","4","5"),nrow=1)
newData=permuteData(dataLine=dataLine,fromCol=2)
}
\keyword{permutation }
