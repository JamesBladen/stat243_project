\name{sampling}
\alias{sampling}
\title{sampling}
\arguments{
  \item{object}{\code{\linkS4class{Cadapt_reject_sample}}
  object}
}
\description{
  Method to sample from s_x.  The basic algorithm is as
  follows: 1. Determine an interval to sample from using
  the weights of integration of the function on each
  interval, computed in the s_x method. 2.  Use inverse CDF
  method to sample from within the selected interval.
  Return the object with new sample.
}

