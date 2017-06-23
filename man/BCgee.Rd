\name{BCgee}
\alias{BCgee}
\alias{print.BCgee}
\alias{summary.BCgee}
\alias{print.summary.BCgee}

\title{
Function to obtain bias-corrected estimates for regression coefficients of a marginal model estimated with generalized estimating equations.
}
\description{
The function takes as an input an object of class \code{gee} and produces an object of class \code{BCgee} which contains the bias-corrected estimates of regression coefficients and further quantities; see \code{details}. 
}
\usage{
BCgee(fit)
}

\arguments{
\item{fit}{
 A fitted model from \code{gee} package. The object must be of class \code{gee}.
}
}

\details{
The function computes bias-corrected estimates of regression coefficients by using quantities in the supplied object of class \code{gee}. The output is an object of class \code{BCgee} which has the same structure of an object of class \code{gee}. The output of the two classes differ because quantities included in the object of class \code{BCgee} are computed with the bias-corrected estimates. 

Methods print and summary are available for objects of class \code{BCgee}.
}

\references{
  Lunardon, N., Scharfstein, D. (2017).
  Comment on "Small sample GEE estimation of regression parameters for longitudinal data".
  \emph{Statistics in Medicine}, <doi:10.1002/sim.7366>. 
}

\author{Lunardon, N. \email{nicola.lunardon@unimib.it}}
\seealso{
  \code{\link{glm}}, \code{\link[gee:gee]{gee}}
}

\examples{
##Cerebrovascular deficiency example
##see page 153 of Diggle, P., Liang, K.-Y., Zeger, S. (1994). Analysis of longitudinal data. 
data(cereb)

if(require(gee)){
fit <- gee(y ~ Period+Drug, id = id, data = cereb, family = binomial(logit), 
corstr = "exchangeable")
fitbc <- BCgee(fit)
##compare coefficients, standard errors, and Wald statistics
summary(fit)$coefficients
summary(fitbc)$coefficients
##compare residuals
fit$residuals
fitbc$residuals
}

##Seizure example from geepack
##see page 166 of Diggle, P., Liang, K.-Y., Zeger, S. (1994). Analysis of longitudinal data. 
data(seizure)
seiz.l <- reshape(seizure,
          varying=list(c("base","y1", "y2", "y3", "y4")),
          v.names="y", times=0:4, direction="long")
     seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time),]
     seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
     seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)


if(require(gee)){
fit <- gee(y ~ offset(log(t)) + x + trt + x:trt, id = id, data=seiz.l, corstr="exchangeable", 
family=poisson(log))
fitbc <- BCgee(fit)
##compare coefficients, standard errors, and Wald statistics
summary(fit)$coefficients
summary(fitbc)$coefficients
##compare residuals
fit$residuals
fitbc$residuals
}
}

