\name{Armitage}
\alias{Armitage}
\title{
Armitage's trend test for the 2x3 genotype table
}
\description{
Armitage's trend test for the 2x3 genotype table.
}
\usage{
Armitage(genotypeStat)
}
\arguments{
  \item{genotypeStat}{
a matrix of genotype frequency with two rows indicating case and control, three column for genotype frequency
(i.e. AA,AG,GG). 
}
}
\details{
The Cochran-Armitage test for trend, is used in categorical data analysis when the aim is to assess for the 
presence of an association between a variable with two categories and a variable with k categories. It modifies
the Pearson chi-squared test to incorporate a suspected ordering in the effects of the k categories of the second
variable. For example, doses of a treatment can be ordered as 'low','medium',and 'high', and we may suspect that
the treatment benefit cannot become smaller as the dose increases. The trend test is ofen used as a genotype-based
test for case-control genetic association studies.
}
\value{
  \item{statistic }{numeric,the statistic of armitage test for trend.}
  \item{pValue }{numeric,the p value of armitage test for trend.}
}
\references{
statgen.org(2007):A derivation for Armitage's trend test for the 2x3 genotype table.
}
\author{
Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\note{
The matrix of genotypeStat must has two rows and three column, otherwise the return value will be NA.
}
\seealso{
  \code{\link{permuteData}}, 
  \code{\link{permuteMatrix}},   
  \code{\link{OR}},   
  \code{\link{Armitage}},
  \code{\link{genotypeStat}},  
  \code{\link{tradPerm}},
  \code{\link{MCPerm}}
}
\examples{
temp=matrix(c(34,56,39,23,45,67),nrow=2,byrow=TRUE)
Armitage(temp)
}
\keyword{Armitage }
