\name{ars}
\alias{ars}
\title{ars: The adapt_reject function}
\usage{
  ars(n_samples, fx, bounds = c(-Inf, Inf),
    guess_of_mode = -999, ...)
}
\arguments{
  \item{n_samples:}{Number of samples desired from
  distribution}

  \item{fx:}{Function to sample from}

  \item{bounds:}{Bounds of function of interest.  The
  default is an unbounded function}
}
\value{
  a vector containing \code{n} points sampled from the f(x)
  distribution
}
\description{
  This calls the class Cadapt_reject_sample and its
  methods.  It returns a vector of samples generated via
  the Adaptive rejective sampling method.
}

