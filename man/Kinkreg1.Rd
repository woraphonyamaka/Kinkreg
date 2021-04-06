\name{Kinkreg1}
\alias{Kinkreg1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kink regression estimation without exogenous varible
}
\description{
The least sqaure estimation of kink regression model without exogenous varible
}
\usage{
Kinkreg1(y,x,Call="optim")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the vector of dependent variable}
  \item{x}{ the matrix of regime-dependent independent variable}
  \item{Call}{ if Call="optim"     == optimization
               if Call="residuals" == retrun the residul series}
}
\details{
The kink regression model without exogenous variable is the nonlinear model in which the continous function is change at the kink point.
}
\value{
\item{Result}{the parameter estimates, standard error, t-stat and p-value}
\item{RSS}{Residual sum Sqaures}
}
\references{
Hansen, B. E. (2017). Regression kink with an unknown threshold. Journal of Business & Economic Statistics, 35(2), 228-240.

Sriboochitta, S., Yamaka, W., Maneejuk, P., & Pastpipatkul, P. (2017). A generalized information theoretical approach to non-linear time series model. In Robustness in Econometrics (pp. 333-348). Springer, Cham.

Maneejuk, P., Yamaka, W., & Sriboonchitta, S. (2019). Does the Kuznets curve exist in Thailand?
A two decades’ perspective (1993–2015). Annals of Operations Research, 1-32.

Maneejuk, P., & Yamaka, W. (2020). An analysis of the impacts of telecommunications technology
and innovation on economic growth. Telecommunications Policy, 44(10), 102038.

Maneejuk, P., & Yamaka, W. (2021). The Impact of Higher Education on Economic Growth in
ASEAN-5 Countries. Sustainability, 13(2), 520.

}
\author{
Woraphon Yamaka
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#### Example Simulation data

set.seed(111)                 # Set seed for reproducibility
k = 1                         #number of dependent variable
n = 500                       #number of observation
r1=1.5                        # Kink parameter
x = matrix(rnorm(n*k,r1,2),ncol=k)
e = rnorm(n,0,sd=1)
x1 = cbind(neg.part(x-r1),pos.part(x-r1))
y=0.5+(0.5*x1[,1])-(1*x1[,2])+e
plot(x,y, col="red", lwd=2)

Kinkreg1(y,x,Call="optim")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
