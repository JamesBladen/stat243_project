\name{ARSpackage}
\alias{ARSpackage}
\alias{ARSpackage}
\docType{package}
\title{
{ARSpackage:} an Adaptive Rejection Sampler
}
\description{
Final project for Statistics 243, an R package that performs adaptive rejection sampling, first proposed by Gilks and Wild in 1992.
}
\details{
\tabular{ll}{
Package: \tab ARSpackage\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-13\cr
Depends: \tab methods, numDeriv\cr
Collate: \tab 'adapt_reject.r', 'ars_methods.r'
}
}
\author{
J. Bladen, L. Felberg, H.W. Tsao, S. Tu
}
\references{
Gilks, Wild, 1992. \url{http://faculty.chicagobooth.edu/hedibert.lopes/teaching/ccis2010/1992GilksWild.pdf}.
}
\keyword{ package, rejection sampling }
\seealso{
\url{https://bitbucket.org/lfelberg/stat243_final_proj} \cr
\url{https://github.com/paciorek/stat243-fall-2013/tree/master/project}
}

\examples{
	# Testing the normal distribution from -Inf to Inf
	n_samples <- ars( 10000, fx = function(x){(1/sqrt(2*pi)*exp((-(x-0)^2)/2))}, bounds=c(-Inf, Inf) )

	# Testing Gamma(2,1) on interval[0.01,Inf]
	sample<-ars(10000,function(x){1/2*x*exp(-x)},c(0.01, Inf))
}
