\name{lower}
\alias{lower}
\title{Cadapt_reject_sample lower}
\arguments{
  \item{object}{\code{\linkS4class{Cadapt_reject_sample}}
  object}

  \item{x}{\code{vector} that we have evaluated h_x for}

  \item{x_star}{\code{numeric} random number for
  adapt/reject}
}
\value{
  l_x_star \code{numeric} the lower hull for \eqn{x^*}
}
\description{
  This calculates the lower hull for \eqn{x^*} which we
  sample from the sampling method for adapt/reject.
}

