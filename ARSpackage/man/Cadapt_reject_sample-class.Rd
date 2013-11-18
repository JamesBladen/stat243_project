\name{Cadapt_reject_sample-class}
\Rdversion{1.1}
\docType{class}
\alias{Cadapt_reject_sample-class}
\alias{error_check,Cadapt_reject_sample-method}
\alias{gen_x,Cadapt_reject_sample-method}
\alias{initialize,Cadapt_reject_sample-method}
\alias{lower,Cadapt_reject_sample-method}
\alias{sample,Cadapt_reject_sample-method}
\alias{show,Cadapt_reject_sample-method}
\alias{update,Cadapt_reject_sample-method}
\alias{upper,Cadapt_reject_sample-method}

\title{Class \code{"Cadapt_reject_sample"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Cadapt_reject_sample", n, h_x, h_prime)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{n}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{h_x}:}{Object of class \code{"function"} ~~ }
    \item{\code{h_prime}:}{Object of class \code{"function"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{error_check}{\code{signature(object = "Cadapt_reject_sample")}: ... }
    \item{gen_x}{\code{signature(object = "Cadapt_reject_sample")}: ... }
    \item{initialize}{\code{signature(.Object = "Cadapt_reject_sample")}: ... }
    \item{lower}{\code{signature(object = "Cadapt_reject_sample")}: ... }
    \item{sample}{\code{signature(object = "Cadapt_reject_sample")}: ... }
    \item{show}{\code{signature(object = "Cadapt_reject_sample")}: ... }
    \item{update}{\code{signature(object = "Cadapt_reject_sample")}: ... }
    \item{upper}{\code{signature(object = "Cadapt_reject_sample")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("Cadapt_reject_sample")
}
\keyword{classes}
