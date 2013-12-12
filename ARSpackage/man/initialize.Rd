\name{initialize}
\alias{initialize}
\title{Cadapt_reject_sample initialization}
\arguments{
  \item{object}{\code{\linkS4class{Cadapt_reject_sample}}
  object}

  \item{n}{\code{numeric} determining the number of samples
  to obtain}

  \item{f_x}{\code{function} for distribution to sample
  from}

  \item{bounds}{\code{vector} of distribution bounds}

  \item{guess_of_mode}{\code{numeric} optional idea of
  where distribution is located}
}
\description{
  A method to intialize the ARS class for sampling.  Will
  store values input from user and will also initialize
  empty arrays for all other slots.
}

