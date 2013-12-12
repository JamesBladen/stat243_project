\name{s_x}
\alias{s_x}
\title{s(x)}
\arguments{
  \item{object}{\code{\linkS4class{Cadapt_reject_sample}}
  object}
}
\description{
  Cadapt_reject_sample method to normalize the upper bounds
  of log(f(x)).  Multiple objective are performed here. The
  most important being the calculation of the abcissa
  vector Z.  Additionally, the weights and exact values of
  the piecewise integration of each interval and the
  normalization factor for the entire upper bound are
  calculated and the x's, their evaluations and their
  derivatives are sorted by x.
}

