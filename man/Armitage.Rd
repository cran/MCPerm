\name{Armitage}
\alias{Armitage}
\title{
Armitage's trend test for the 2x3 genotype table
}
\description{
Armitage's trend test for the 2x3 genotype table.
}
\usage{
Armitage(case_11, case_12, case_22, control_11, control_12, control_22)
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
}
\details{
The Cochran-Armitage test for trend, is used in categorical data analysis when the aim is to assess
for the presence of an association between a variable with two categories and a variable with k
categories. It modifies the Pearson chi-squared test to incorporate a suspected ordering in the effects
of the k categories of the second variable. The trend test is ofen used as a genotype-based test for 
case-control genetic association studies.
}
\value{
  \item{statistic }{numeric, the statistic of armitage test for trend.}
  \item{pValue }{numeric, the p value of armitage test for trend.}
}
\references{
Armitage, P(1955): Tests for Linear Trends in Proportions and Frequencies.

statgen.org(2007): A derivation for Armitage's trend test for the 2x3 genotype table.
}
\author{
Lanying Zhang and Yongshuai Jiang <jiangyongshuai@gmail.com>
}
\seealso{
    \code{\link{Armitage.TradPerm}},
    \code{\link{Armitage.MCPerm}},
    \code{\link{OR}},
    \code{\link{OR.TradPerm}},
	\code{\link{OR.MCPerm}},
	\code{\link{permuteData}},
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
case_11=4
case_12=1
case_22=1
control_11=3
control_12=0
control_22=0
Armitage(case_11,case_12,case_22,control_11,control_12,control_22)
}
\keyword{ Armitage }
