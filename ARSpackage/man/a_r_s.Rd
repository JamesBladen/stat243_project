\name{a_r_s}
\alias{a_r_s}
\title{The adapt_reject function}
\usage{
  a_r_s(n_samples, log_fx, log_fx_prime, ...)
}
\arguments{
  \item{n_samples}{Number of samples desired from
  distribution}

  \item{log_fx}{Log of function to sample from}

  \item{log_fx_prime}{First derivative of log of function
  to sample from}
}
\value{
  S4 \code{adapt_reject_sample} object; a vector containing
  \code{n} points sampled from the f(x) distribution
}
\description{
  This calls the class Cadapt_reject_sample and its
  methods.
}

