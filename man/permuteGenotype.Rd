\name{permuteGenotype}
\alias{permuteGenotype}
\title{
Fill the numerics of 2*3 table when fixed the row and column totals
}
\description{
Fill the numerics of 2*3 table when fixed the row and column totals.
}
\usage{
permuteGenotype(case_11, case_12, case_22, control_11, control_12, control_22, n)
}
\arguments{
  \item{case_11}{
a non-negative integer, the frequency of genotype "allele1/allele1" in case samples. 
}
  \item{case_12}{
a non-negative integer, the frequency of genotype "allele1/allele2" in case samples. 
}
  \item{case_22}{
a non-negative integer, the frequency of genotype "allele2/allele2" in case samples. 
}
  \item{control_11}{
a non-negative integer, the frequency of genotype "allele1/allele1" in control samples. 
}
  \item{control_12}{
a non-negative integer, the frequency of genotype "allele1/allele2" in control samples. 
}
  \item{control_22}{
a non-negative integer, the frequency of genotype "allele2/allele2" in control samples. 
}
  \item{n}{
an integer specifying the number to generate.
}
}
\details{
Fill the numerics of 2*3 table when fixed the row and column totals, by twice generating 
random numbers for the hypergeometric distribution, based on the putting data.
}
\value{
  \item{perm_case_11    }{a vector with 'n' elements.}
  \item{perm_case_12    }{a vector with 'n' elements.}
  \item{perm_case_22    }{a vector with 'n' elements.}
  \item{perm_control_11 }{a vector with 'n' elements.}
  \item{perm_control_12 }{a vector with 'n' elements.}
  \item{perm_control_22 }{a vector with 'n' elements.}
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
	\code{\link{permuteData}},
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
permuteGenotype(case_11=23,case_12=0,case_22=34,control_11=0,control_12=34,control_22=45,n=5)
}
\keyword{ permutation }
\keyword{ genotype }
