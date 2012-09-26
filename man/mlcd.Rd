\name{mlcd}
\alias{mlcd}
\docType{data}
\title{
Multivariate Longitudinal Count Data (MLCD)
}
\description{
A data frame with 2000 observations on the following 6 variables. MLCD is a simulated bivariate longitudinal count dataset assuming there are 500 subjects in the study whose data are collected at 4 equally-spaced time points. 
}
\usage{data(mlcd)}
\format{
  A data frame with 2000 observations on the following 6 variables.
  \describe{
    \item{\code{ID}}{a numeric vector for subject ID}
    \item{\code{resp1}}{a numeric vector for the first longitudinal count response}
    \item{\code{resp2}}{a numeric vector for the second longitudinal count response}
    \item{\code{X}}{a numeric vector for the covariate, X}
    \item{\code{time}}{a numeric vector for the time point at which observations are collected}
    \item{\code{X.time}}{a numeric vector for the interaction between X and time}
  }
}
\details{
The covariates, X and time are the standardized values indeed. The related interaction is calculated by using these standardized values. X is a time-independent covariate.
}
%\source{
%
%}
\references{
Asar, O. (2012). \emph{On multivariate longitudinal binary data models and their applications in forecasting}. MS Thesis, Middle East Technical University

Erhardt, V. (2009). corcounts: Generate Correlated Count Random Variable. R package version 1.4. 
URL http://CRAN.R-project.org/package=corcounts.
}
\examples{
data(mlcd)
plot(mlcd$X,mlcd$resp1)
}
\keyword{datasets}
