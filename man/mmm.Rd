\name{mmm}
\alias{mmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mmm: An R Package for Multivariate marginal models
}
\description{
fits multivariate marginal models for both continous and discrete responses
}
\usage{
mmm(formula, id, data = NULL, R = NULL, b = NULL, tol = 0.001, maxiter = 25, 
family = "gaussian", corstr = "independence", Mv = 1, silent = TRUE, 
scale.fix = FALSE, scale.value = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula expression, see the examples given below.
}
  \item{id}{
a vector for identification of the clusters.
}
  \item{data}{
an optional data frame.
}
  \item{R}{
a user specified square matrix for the working correlation matrix, appropriate when \code{corstr="fixed"}.
}
  \item{b}{
user specified initials for the parameter estimates.
}
  \item{tol}{
the tolerance which specifies the convergency of the algorithm. 
}
  \item{maxiter}{
the maximum number of iterations to be consumed by the algorithm.
}
  \item{family}{
an object which defines the link and variance function. The possible choices are same with the ones in the \code{"gee"} package. For details see the \code{\link{gee}} documentation. Note that \code{family=binomial} handles multivariate longitudinal binary data, \code{family=poisson} handles multivariate longitudinal count data, \code{family=gaussian} handles multivariate longitudinal (normal type) continous data and \code{family=gamma} handles multivariate longitudinal (gamma type) continous data. 
}
  \item{corstr}{
a character string which defines the structure of the working correlation matrix. For details see the \code{\link{gee}} documentation.
}
  \item{Mv}{
specifies the lag value, e.g. specification of \code{"corstr=AR-M"} and \code{"Mv=1"} indicates AR(1).
}
  \item{silent}{
a logial variable which decides the print of the iterations.
}
  \item{scale.fix}{
a logical variable for fixing the scale parameter to a user specified value.
}
  \item{scale.value}{
a user specified scale parameter value, appropriate when \code{scale.fix=TRUE}.
}
}

\value{
Returns an onject of the results. See the examples given below.
}
\references{
Asar, O. (2012). \emph{On multivariate longitudinal binary data models and their applications in forecasting}. MS Thesis, Middle East Technical University

Liang, K. L., Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika}, \bold{73}, 13-22.

Shelton, B. J., Gilbert, G. H., Liu, B., Fisher, M. (2004). A SAS macro for the analysis of multivariate longitudinal binary outcomes. \emph{Computer Methods and Programs in Biomedicine}, \bold{76}, 163-175. 

Zeger, S. L., Liang, K. L (1986). Longitudinal data analysis for discrete and continous outcomes. \emph{Biometrics}, \bold{42}, 121-130.
}
\author{
Ozgur Asar, Ozlem Ilk
}
\note{
Version 1.1.
}

\seealso{
\code{\link{gee}}
}
\examples{
## Binary data example
data(mscm)
fit1<-mmm(formula=cbind(mscm$stress,mscm$illness)~mscm$married+mscm$education+
mscm$employed+mscm$chlth+mscm$mhlth+mscm$race+mscm$csex+mscm$housize+
mscm$bstress+mscm$billness+mscm$week,id=mscm$id,family=binomial,corstr="exchangeable")
summary(fit1)

## Count data example

data(mlcd)
fit2<-mmm(formula=cbind(mlcd$resp1,mlcd$resp2)~mlcd$X+mlcd$time+
mlcd$X.time,id=mlcd$ID,family=poisson,corstr="unstructured")
summary(fit2)

## Continuous data example

data(mlgd)
fit3<-mmm(formula=cbind(mlgd$resp1,mlgd$resp2)~mlgd$X+mlgd$time+
mlgd$X.time,id=mlgd$ID,family=gaussian,corstr="unstructured")
summary(fit3)
}

\keyword{ generalized estimating equations }

