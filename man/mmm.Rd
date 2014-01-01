\name{mmm}
\alias{mmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
an R function to fit the multivariate marginal models to analyze multivariate longitudinal data
}
\description{
fits multivariate marginal models to analyze multivariate longitudinal data, for both continuous and discrete responses
}
\usage{
mmm(formula, id, data = NULL, correlation = NULL, initEstim = NULL, tol = 0.001, 
maxiter = 25, family = "gaussian", corStruct = "independence", Mv = 1, 
silent = TRUE, scale.fix = FALSE, scale.value = 1)
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
  \item{correlation}{
a user specified square matrix for the working correlation matrix, appropriate when \code{corstr="fixed"}.
}
  \item{initEstim}{
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
  \item{corStruct}{
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
Asar, O. (2012). \emph{On multivariate longitudinal binary data models and their applications in forecasting}. MS Thesis, Middle East Technical University.
Available online at http://www.lancaster.ac.uk/pg/asar/thesis-Ozgur-Asar

Asar, O., Ilk, O. (2013). mmm: an R package for analyzing multivariate longitudinal data with multivariate marginal models, 
\emph{Computer Methods and Programs in Biomedicine}, \bold{112}, 649-654.

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
#########################
## Binary data example ##
#########################

data(motherStress)
fit1<-mmm(formula=cbind(stress,illness)~married+education+
employed+chlth+mhlth+race+csex+housize+bstress+billness+
week,id=motherStress$id,data=motherStress,family=binomial,corStruct="exchangeable")
summary(fit1)

########################
## Count data example ##
########################

## First we illustrate how the data set is simulated
## Then the R script to analyze the data set by mmm is given
## Note: no need to run the script to generate the data set, unless of interest

\dontrun{
### Generating the data by the help of 'corcounts' package

# loading the package 'corcounts'
library("corcounts")

# setting the seed to 12
set.seed(12)

# number of subjects in the study
n1 <- 500

# defining the response and covariate families (Poi indicates Poisson distribution)
margins <- c("Poi","Poi","Poi","Poi","Poi","Poi","Poi","Poi","Poi")

# the means of the responses and covariate. while 5 and 8 are the means of the responses
# 20 is the mean of the time independent covariate
mu <- c(5, 8, 20, 5, 8, 5, 8, 5, 8)

# the correlation structure which 'corcounts' uses to generate correlated count data
# (unstr indicates unstructured correlation structure)
corstr <- "unstr"

# the correlation matrix corcounts assumes the correlated count data have
corpar<-matrix(c(1,0.4,0.6,0.9,0.37,0.8,0.34,0.7,0.31,
0.4,1,0.6,0.37,0.9,0.34,0.8,0.31,0.7,
0.6,0.6,1,0.6,0.6,0.6,0.6,0.6,0.6,
0.9,0.37,0.6,1,0.4,0.9,0.37,0.8,0.34,
0.37,0.9,0.6,0.4,1,0.37,0.9,0.34,0.8,
0.8,0.34,0.6,0.9,0.37,1,0.4,0.9,0.37,
0.34,0.8,0.6,0.37,0.9,0.4,1,0.37,0.9,
0.7,0.31,0.6,0.8,0.34,0.9,0.37,1,0.4,
0.31,0.7,0.6,0.34,0.8,0.37,0.9,0.4,1),ncol=9,byrow=T)

# generating the correlated count data by 'rcounts' function avaiable in 'corcounts'
data1 <- rcounts(N=n1,margins=margins,mu=mu,corstr=corstr,corpar=corpar)

### The reconstruction of the generated correlated count data to 
### the longitudinal data (long) format

# seperating the bivariate responses measured at the first time 
# point and the time independent covariate 
time11<-data1[,1:3]

# seperating the bivariate responses measured at the second time
# point and combining them with the time independent covariate
time12<-cbind(data1[,4:5],data1[,3])

# seperating the bivariate responses measured at the third time
# point and combining them with the time independent covariate
time13<-cbind(data1[,6:7],data1[,3])

# seperating the bivariate responses measured at the fourth time
# point and combining them with the time independent covariate
time14<-cbind(data1[,8:9],data1[,3])

# combining the data for all the time points
data12<-rbind(time11,time12,time13,time14)

# constructing the time variable
time1<-matrix(rep(seq(1:4),each=n1))

# constructing the id variable 
id1<-matrix(rep(seq(1:n1),4))

# combining the id of the subjects, the simulated data and the time variable
data13<-cbind(id1,data12,time1)

# reconstructing the data subject by subject which 'mmm' expects it has 
data14<-NULL
for (i in 1:n1) data14<-rbind(data14,data13[data13[,1]==i,])

### Data manipulations on the covariates

# taking natural logarithm of the time independent covariate
data14[,4]<-log(data14[,4])

# standardizing time variable
data14[,5]<-scale(data14[,5])

# adding the interaction of the time independent covariate 
# and time as a new covariate
multiLongCount<-as.data.frame(cbind(data14,data14[,4]*data14[,5]))
names(multiLongCount)<-c("ID","resp1","resp2","X","time","X.time")
}

### R script to analyze the count data set
### It is already stored in mmm pacakge

data(multiLongCount)
fit2<-mmm(formula=cbind(resp1,resp2)~X+time+X.time,
id=multiLongCount$ID,data=multiLongCount,family=poisson,corStruct="unstructured")
summary(fit2)

#############################
## Continuous data example ##
#############################

## First we illustrate how the data set is simulated
## Then the R script to analyze the data set by mmm is given
## Note: no need to run the script to generate the data set, unless of interest

\dontrun{
### Generating the data by the help of mvtnorm package

# loading package 'mvtnorm' 
library("mvtnorm")

# number of subjects in the study
n2<-500

# setting the seed to 12
set.seed(12)

# specifying the correlation matrix which 'mvtnorm' assumes the correlated data have
cormat<-matrix(c(1,0.4,0.6,0.9,0.37,0.8,0.34,0.7,0.31,
0.4,1,0.6,0.37,0.9,0.34,0.8,0.31,0.7,
0.6,0.6,1,0.6,0.6,0.6,0.6,0.6,0.6,
0.9,0.37,0.6,1,0.4,0.9,0.37,0.8,0.34,
0.37,0.9,0.6,0.4,1,0.37,0.9,0.34,0.8,
0.8,0.34,0.6,0.9,0.37,1,0.4,0.9,0.37,
0.34,0.8,0.6,0.37,0.9,0.4,1,0.37,0.9,
0.7,0.31,0.6,0.8,0.34,0.9,0.37,1,0.4,
0.31,0.7,0.6,0.34,0.8,0.37,0.9,0.4,1),ncol=9,byrow=T)

# variances of the responses and time independent covariate
# while 0.97 and 1.1 correspond to the variances of the bivariate responses
# 4 corresponds to the variance of the time independent covariate
variance<-c(0.97,1.1,4,0.97,1.1,0.97,1.1,0.97,1.1)

# constructing the (diaonal) standard deviation matrix 
std<-diag(sqrt(variance),9)

# constructing the variance covariance matrix, sigma
sigma<-std%*%cormat%*%std

# generating the correlated continuous data utilizing 'rmvnorm' function available
# in 'mvtnorm'; method="svd" indicates use of 'singular value decomposition method
data2<-rmvnorm(n2,mean = rep(0,nrow(sigma)),sigma=sigma,method="svd")


### The reconstruction of the generated correlated continuous data to the 
### longitudinal data (long) format

# seperating the bivariate responses measured at first time point 
# and the time independent covariate
time21<-data2[,1:3]

# seperating the bivariate responses measured at second time point
# and combining them with the time independent covariate
time22<-cbind(data2[,4:5],data2[,3])

# seperating the bivariate responses measured at third time point
# and combining them with the time independent covariate
time23<-cbind(data2[,6:7],data2[,3])

# seperating the bivariate responses measured at fourth time point
# and combining them with the time independent covariate
time24<-cbind(data2[,8:9],data2[,3])

# combining the data for all the time points
data22<-rbind(time21,time22,time23,time24)

# constructing the time variable
time2<-matrix(rep(seq(1:4),each=n2))

# constructing the id variable
id2<-matrix(rep(seq(1:n2),4))

# combining the id of the subjects, the generated data and the time variable
data23<-cbind(id2,data22,time2)

# reconstructing the data subject by subject which 'mmm' expects it has
data24<-NULL
for (i in 1:n2) data24<-rbind(data24,data23[data23[,1]==i,])

### Data manipulations on the covariates

# standardizing the time variable
data24[,5]<-scale(data24[,5])

# adding the interaction of the time independent covariate 
# and time as a new covariate
multiLongGaussian<-as.data.frame(cbind(data24,data24[,4]*data24[,5]))
names(multiLongGaussian)<-c("ID","resp1","resp2","X","time","X.time")
}


### R script to analyze the continuous data set
### It is already stored in mmm pacakge

data(multiLongGaussian)
fit3<-mmm(formula=cbind(resp1,resp2)~X+time+X.time,
id=multiLongGaussian$ID,data=multiLongGaussian,family=gaussian,corStruct="unstructured")
summary(fit3)
}

\keyword{ generalized estimating equations }

