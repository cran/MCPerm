\name{exampleData}
\alias{exampleData}
\docType{data}
\title{
Data Example for Function tradPerm and genotypeStat
}
\description{
Data Example for Function tradPerm and genotypeStat.
}
\usage{data(exampleData)}
\format{
  A data frame with 6 observations on the following 48 variables.
  \describe{
    The data frame with 48 Column:
	  the first column: snp id;
      47 column: genotype data,	e.g. A_A,A_G,G_G,?_?.  
  }
}
\details{
 Data Example for Function tradPerm and genotypeStat.
}
\references{
Edgington. E.S. (1995) Randomization tests, 3rd ed. New York: Marcel-Dekker
}
\examples{
## import example data
data(exampleData)
data1=exampleData[1,,drop=FALSE]
data2=exampleData[2,,drop=FALSE]
result1=tradPerm(data2,data1,fromCol=2,naString="?_?",sep="_",testMethod="fisher",repeatNum=10)
result2=genotypeStat(data2,data1,fromCol=2,naString="?_?",sep="_")
}
\keyword{datasets}
