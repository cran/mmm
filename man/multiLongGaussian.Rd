\name{multiLongGaussian}
\alias{multiLongGaussian}
\docType{data}
\title{
Multivariate Longitudinal Continuous (Gaussian) Data
}
\description{
A data frame with 2000 observations on the following 6 variables. multiLongGaussian is a simulated bivariate longitudinal continuous dataset assuming there are 500 subjects in the study whose data are collected at 4 equally-spaced time points. 
}
\usage{data(multiLongGaussian)}
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
The covariates, X and time are the standardized values indeed. The related interaction is calculated by using these standardized values. X is a time-independent covariate. The 
R script to generate the data set is given in the Examples section of the mmm function.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Asar, O. (2012). \emph{On multivariate longitudinal binary data models and their applications in forecasting}. MS Thesis, Middle East Technical University

Genz, A., Bretz, F., Miwa, T., Mi, X., Leisch, F., Scheipl, F., Hothorn, T. (2011). mvtnorm: Multivariate Normal and t Distributions. R package version 0.9-96. 
URL http://CRAN.R-project.org/package=mvtnorm.
}
\examples{
data(multiLongGaussian)
plot(multiLongGaussian$X,multiLongGaussian$resp1)
}
\keyword{datasets}
