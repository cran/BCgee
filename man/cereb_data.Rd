\name{cereb}
\alias{cereb}
\non_function{}
\title{Cerebrovascular Deficiency}
\usage{data(cereb)}
\description{
The \code{cereb} data frame has 134 rows and 4 columns. The dataset consists of safety data from a crossover trial on the disease cerebrovascular deficiency. In this two-period crossover trial, comparing the effects of active drug to placebo, 67 patients were randomly allocated to the two treatment sequences, with 34 patients receiving placebo followed by active treatment, and 33 patients receiving active treatment followed by placebo. The response variable is binary, indicating whether an electrocardiogram (ECG) was abnormal (Y=1) or normal (Y=0). Each patient has a bivariate binary response vector.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{patient's unique id number}
    \item{Period}{Period in the trial, 0=Period 1, 1=Period 2}
    \item{Drug}{Treatment, 0=Placebo, 1=Active Drug}
    \item{y}{ECG Response, 0=Normal, 1=Abnormal}    
  }
}
\references{
  Jones, B. and Kenward, M.G. (1989). Design and Analysis of Cross-over Trials. London: Chapman and Hall/CRC Press.

  Diggle, P.J., Liang, K.Y., and Zeger, S.L. (1994). Analysis of Longitudinal Data. Clarendon Press.
}
\examples{
##see help page of function BCgee
}
\keyword{datasets}
