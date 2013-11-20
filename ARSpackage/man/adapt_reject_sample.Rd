\name{Cadapt_reject_sample}
\alias{Cadapt_reject_sample}
\title{The adapt_reject class}
\value{
  S4 \code{adapt_reject_sample} object; a vector containing
  \item{n} points sampled from the f(x) distribution
}
\description{
  This file consists of the adapt_reject_sample class, the
  function that calls it and the functions it uses.
}
\details{
  It performs an adaptive rejection sampling process as
  proposed by Wild and Gilks in 1992.
}
\note{
  1. Initialize i) x1, x2 ii) inputs: h(x) and h'(x), n
  (number of points to sample), optional: domain etc iii)
  error checks: make sure that the function is concave up
  and the function lies within U(x) and L(x).  Check that
  x1 has a positive slope and X2 has a negative slope.
  Check that the sample size is positive and an integer. 2)
  Objects/methods: i) U(x) and S(x): z(x), equations for
  tangent lines ii) List of x points iii) list of sampled
  points iv) l(x) v) sample function from s(x) and uniform
  random number vi) update steps vii) error checking

  Current questions: 1.  How do we draw a random number
  from sk(x) i) calculate the area under each piece (Sk(x))
  ii) divide by total area (Stot(x)) iii) weights <-
  Sk(x)/Stot(x) iv) sample(1:k+1 with weights) -> select
  piece v) rejection sample within the piece

  OR maybe a package?? spatstat with rpoint

  2.  How do we find initial points for x1 and x2?? All we
  know right now is that they need to encompass the max?
  One needs pos deriv and one needs neg method 1: gen
  random number and calculate h_prime method 2: find 2
  stdevs from mean, check if they fit criteria
}
\section{Slots}{
  \describe{ \item{\code{n}:}{Variable of class
  \code{"numeric"}, n, containing the number of points to
  sample} \item{\code{h_x}:}{Function of class
  \code{"function"}, containing the log(f(x)) to sample
  from.} \item{\code{h_prime}:}{Function of class
  \code{"function"}, containing the first derivative
  log(f(x)) to sample from.} \item{\code{x}:}{Variable of
  class \code{"vector"}, containing points used to draw
  lines.} \item{\code{z}:}{Variable of class
  \code{"vector"}, containing abscissae of upper bound
  function.} \item{\code{output}:}{Variable of class
  \code{"vector"}, containing sampled points to return to
  user.} }
}
\author{
  J. Bladen, L. Felberg, H.W. Tsao, S. Tu
}
\references{
  \url{http://faculty.chicagobooth.edu/hedibert.lopes/teaching/ccis2010/1992GilksWild.pdf}
}

