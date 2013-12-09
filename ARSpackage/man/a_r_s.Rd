\name{a_r_s}
\alias{a_r_s}
\title{The adapt_reject function}
\usage{
  a_r_s(n_samples, fx, bounds = c(-Inf, Inf), ...)
}
\arguments{
  \item{n_samples:}{Number of samples desired from
  distribution}

  \item{fx:}{Function to sample from}

  \item{bounds:}{Bounds of function of interest.  The
  default is an unbounded function}
}
\value{
  S4 \code{adapt_reject_sample} object; a vector containing
  \code{n} points sampled from the f(x) distribution
}
\description{
  This calls the class Cadapt_reject_sample and its
  methods.  The vector of samples is accessible via
  \var{ans output}.
}

