\name{gen_x}
\alias{gen_x}
\title{gen_x}
\arguments{
  \item{object}{\code{\linkS4class{Cadapt_reject_sample}}
  object}
}
\description{
  Cadapt_reject_sample method for generating first two
  points.  If the distribution is unbounded, then find the
  function's mode and pick points surrounding it.  If it's
  bounded on one side, we use the bound given and search
  until we find a point that corresponds to the opposite
  end of the domain with respect to their derivatives.  If
  bounded on both sides, use given bounds.
}

