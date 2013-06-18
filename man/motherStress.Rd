\name{motherStress}
\alias{motherStress}
\docType{data}
\title{
Mother's Stress and Children's Morbidity Study
}
\description{
 A data frame with 2004 observations on the following 14 variables. motherStress is a longitudinal dataset which includes daily information of the participants. There are 167 mothers and children enrolled in the study. 
}
\usage{data(motherStress)}
\format{
The details of the columns of the data frame are given below.
  \describe{
    \item{\code{id}}{a vector for subject id}
    \item{\code{stress}}{a vector for mother's stress at time t:1=presence, 0=absence}
    \item{\code{illness}}{a vector for children's illness at time t: 1=presence, 0=absence}
    \item{\code{married}}{a vector for marriage status of mother: 1=married, 0=other}
    \item{\code{education}}{a vector for mother's education level: 0=high school or less, 1=high school graduate}
    \item{\code{employed}}{a numeric vector for mother's employment status: 1=employed, 0=unemployed}
    \item{\code{chlth}}{a vector for children's health status at baseline: 0=very poor/poor, 1=fair, 2=good, 3=very good}
    \item{\code{mhlth}}{a vector for mother's health status at baseline: 0=very poor/poor, 1=fair, 2=good, 3=very good}
    \item{\code{race}}{a vector for child's race: 1=non-white, 0=white}
    \item{\code{csex}}{a vector for child's gender: 1=female, 0=male}
    \item{\code{housize}}{a vector for the size of the household: 0=2-3 people, 1=more than 3 people}
    \item{\code{bstress}}{a vector for the baseline stress for the period of day 1 to 16; calculated as the mean of the stress status of the subjects in the period of day 1 to 16}
    \item{\code{billness}}{a vector for the baseline illness for the period of day 1 to 16; calculated as the mean of the illness status of the subjects in the period of day 1 to 16}
    \item{\code{week}}{a numeric vector for time: (day-22)/7}
  }
}
\details{
The original data contains the information of the mothers and children in the study for 28 days. Because of the weak serial correlation in the period of day 1 to 16, it is ignored. Only the period of day 17 to 28 is included here. To catch the specific characteristic of the mothers and children, the averages of the stress and illness status of them are added as new covariates; bstress and billness. While the covariates have no missing observation, responses have very low percentages of missing values, 0.97% and 1.42% for stress and illness, respectively. Since the scope of \code{mmm} is beyond missing data, those missing values are imputed by occasion mode imputation.
}
\source{
http://faculty.washington.edu/heagerty/Books/AnalysisLongitudinal/datasets.html
}
\references{
Alexander, C. S., Markowitz, R. (1986). Maternal employment and use of pediatric clinic services. \emph{Medical Care}, \bold{24(2)}, 134-147.

Diggle, P. J., Heagerty, P., Liang, K. Y., Zeger, S. L. (2002). \emph{Analysis of Longitudinal Data}. New York: Oxford University Press.

Zeger, S. L., Liang, K. L (1986). Longitudinal data analysis for discrete and continous outcomes. \emph{Biometrics}, \bold{42}, 121-130.

}
\examples{
data(motherStress)
head(motherStress,10)
require(graphics)
mosaicplot(~motherStress$employed+motherStress$housize+motherStress$stress,color=TRUE)
}
\keyword{datasets}

