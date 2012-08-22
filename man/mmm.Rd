\name{mmm}
\alias{mmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit multivariate marginal models with response specific parameter
}
\description{
fits multivariate marginal models with response specific parameter for both continous and discrete responses and produces both the multivariate and univariate marginal modeling results
}
\usage{
mmm(data, nresp, family = "gaussian", corstr = "independence",
coefnames = NULL, tol = 0.001, maxiter = 25, Mv = 1, 
silent = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame. Id of the subjects should be placed in the first column, multivariate responses should be placed in the columns next to id variable and the covariates which are to be included in the model should be placed in the columns next to the multivariate responses.
}
  \item{nresp}{
number of multivariate responses.
}
  \item{family}{
an object which defines the link and variance function. The possible choices are same with the ones in the \code{"gee"} package. For details see the \code{\link{gee}} documentation. Note that \code{family=binomial} handles multivariate longitudinal binary data, \code{family=poisson} handles multivariate longitudinal count data, \code{family=gaussian} handles multivariate longitudinal (normal type) continous data and \code{family=gamma} handles multivariate longitudinal (gamma type) continous data. 
}
  \item{corstr}{
a character string which defines the structure of the working correlation matrix. For details see the \code{\link{gee}} documentation.
}
  \item{coefnames}{
a list of names of the coefficients which provides better output. 
}
  \item{tol}{
the tolerance which specifies the convergency of the algorithm. 
}
  \item{maxiter}{
the maximum number of iterations to be consumed by the algorithm.
}
  \item{Mv}{
specifies the lag value, e.g. specification of \code{"corstr=AR-M"} and \code{"Mv=1"} indicates AR(1).
}
  \item{silent}{
a logial variable which decides the print of the iterations.
}
}
\details{
The \code{mmm} package utilizes the \code{\link{gee}} package within. Therefore, installation and loading of that package is mandatory. \code{mmm} loads it automatically. Note that the \code{mmm} produces the univariate marginal modeling results as well as the multivariate results.}
\value{
\item{list1 }{A list which includes the output of the multivariate marginal modeling and univariate marginal modeling at the same time. For the details of manipulation of the output, see the example below.}
}
\references{
Asar, O. (2012). \emph{On multivariate longitudinal binary data models and their applications in forecasting}. MS Thesis, Middle East Technical University.

Liang, K. L., Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika}, \bold{73}, 13-22.

Shelton, B. J., Gilbert, G. H., Liu, B., Fisher, M. (2004). A SAS macro for the analysis of multivariate longitudinal binary outcomes. \emph{Computer Methods and Programs in Biomedicine}, \bold{76}, 163-175. 

Zeger, S. L., Liang, K. L (1986). Longitudinal data analysis for discrete and continous outcomes. \emph{Biometrics}, \bold{42}, 121-130.

}
\author{
Ozgur Asar, Ozlem Ilk
}
\note{
This is the version 1.0 of this user documentation file.
}

\seealso{
\code{\link{gee}}
}
\examples{
data(mscm.mmm)
mscm<-mscm.mmm
coefnames<-c("str.intercept","str.married","str.education",
"str.employed","str.chlth","str.mhlth","str.race",
"str.csex","str.housize","str.bstress","str.billness",
"str.week","ill.intercept","ill.married","ill.education",
"ill.employed","ill.chlth","ill.mhlth","ill.race",
"ill.csex","ill.housize","ill.bstress","ill.billness",
"ill.week")
mmm1<-mmm(data=mscm,nresp=2,family=binomial,corstr='exchangeable',
coefnames=coefnames)
#coefficients and the related statistics from the 
#multivariate output
round(((mmm1$multiv)$multivout)$coef,4)
#coefficients and the related statistics from the 
#univariate output, for the first response (stress)
round((((mmm1$univ)$univout)[[1]])$coef,4)
#coefficients and the related statistics from the 
#univariate output, for the second response (illness)
round((((mmm1$univ)$univout)[[2]])$coef,4)
}
\keyword{ package}
\keyword{ generalized estimating equations}
