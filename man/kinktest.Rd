\name{kinktest}
\alias{kinktest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Testing the kink effect
%%  ~~function to do ... ~~
}
\description{
Hypothesis testing for kink point in kink regression. The null hyphothesis is  linear regression against the alternative hyphothesis of kink regression.
}
\usage{
kinktest(y,x,level=0.90, boot=10,search=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the vector of dependent variable}
  \item{x}{ the matrix of regime-dependent independent variable}
  \item{level}{ Confidence interval}
  \item{boot}{ number of boostrapping}
   \item{search}{ range of grid search}
}
}
\details{
As shown by Hansen (1996), it is simple
to simulate approximations to (10) using a multiplier bootstrap, and thus asymptotically valid p-values can be calculated. The following is his recommended algorithm. (Theorem 3 of Hansen (1996) shows that the algorithm produces asymptotically first-order correct p-values under the conditions of Theorem 1.)
}
\value{

\item{Wald.test }{Wald statistic}
\item{p-value }{p-value}

}
\references{
Hansen, B. E. (1996). Inference when a nuisance parameter is not identified under the null hypothesis. Econometrica: Journal of the econometric society, 413-430.
Hansen, B. E. (2017). Regression kink with an unknown threshold. Journal of Business & Economic Statistics, 35(2), 228-240.
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
set.seed(111)  # Set seed for reproducibility
k = 1                         #number of dependent variable
n = 500                       #number of observation
r1=1.5                        # Kink parameter
x = rnorm(n,r1,sd=1)
e = rnorm(n,0,sd=1)
x1 = cbind(neg.part(x-r1),pos.part(x-r1))
y=0.5+(0.5*x1[,1])-(1*x1[,2])+e

# Ho: Linear
# Ha: Kink
kinktest(y,x,level=0.90, boot=10,search=0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
