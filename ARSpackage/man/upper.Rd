\name{upper}
\alias{upper}
\title{Cadapt_reject_sample upper}
\arguments{
  \item{object}{\code{\linkS4class{Cadapt_reject_sample}}
  object}

  \item{x}{\code{vector} that we have evaluated h_x for}

  \item{x_star}{\code{numeric} random number for
  adapt/reject}

  \item{z}{\code{vector} abscissa of all points}
}
\value{
  u_x_star \code{numeric} the upper hull for \eqn{x^*}
}
\description{
  This calculates the upper hull for \eqn{x^*} which we
  sample from the sampling method for adapt/reject.
}

