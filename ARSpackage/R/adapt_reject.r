######################################
######################################

#' The adapt_reject function
#'
#' This calls the class Cadapt_reject_sample and its methods.
#'
#' @param n_samples Number of samples desired from distribution
#' @param log_fx Log of function to sample from
#' @param log_fx_prime First derivative of log of function to sample from
#' @return S4 \code{adapt_reject_sample} object; a vector containing \code{n} points sampled from the f(x) distribution

a_r_s <- function( n_samples, log_fx, bounds=c(-Inf, Inf), ... ){
  
  # Initialize new ARS class
  ars_class <- new( "Cadapt_reject_sample", n=n_samples, h_x = log_fx, bounds, ... )
  ars_class <- gen_x( ars_class )
  
  print( ars_class@output )
  
  # While we do not have enough n samples, continue to sample
  while( length( ars_class@output < n_samples ) ){
    
    ars_class <- s_x( ars_class )
    ars_class <- sample( ars_class )
    ars_class <- update( ars_class )
    
  }
  
  return( ars_class )
}

######################################
######################################
#' The adapt_reject class
#'
#' This class contains all the methods used to perform an AR sampling.  
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{n}:}{Variable of class \code{"numeric"}, n, containing the number of points to sample}
#'    \item{\code{h_x}:}{Function of class \code{"function"}, containing the log(f(x)) to sample from.}
#'    \item{\code{n}:}{Variable of class \code{"numeric"}, n, containing the bounds of the function}
#'    \item{\code{x}:}{Variable of class \code{"vector"}, containing points used to draw lines.}
#'    \item{\code{z}:}{Variable of class \code{"vector"}, containing abscissae of upper bound function.}
#'    \item{\code{output}:}{Variable of class \code{"vector"}, containing sampled points to return to user.}
#'  }
#'
#' @note  1. Initialize
#'   i) x1, x2
#'   ii) inputs: h(x) and h'(x), n (number of points to sample), optional: domain etc
#'   iii) error checks: make sure that the function is concave up and the function lies within U(x) and L(x).  
#'         Check that x1 has a positive slope and X2 has a negative slope. Check that the sample size is positive and an integer.
#' 2) Objects/methods:
#'   i) U(x) and S(x): z(x), equations for tangent lines
#'   ii) List of x points
#'   iii) list of sampled points
#'   iv) l(x)
#'   v) sample function from s(x) and uniform random number
#'   vi) update steps
#'   vii) error checking
#'   
#' Current questions:
#'  1.  How do we draw a random number from sk(x)
#'          i) calculate the area under each piece (Sk(x))
#'          ii) divide by total area (Stot(x))
#'          iii) weights <- Sk(x)/Stot(x)
#'          iv) sample(1:k+1 with weights) -> select piece
#'          v) rejection sample within the piece
#'          
#'          OR maybe a package?? spatstat with rpoint
#'          
#'  2.  How do we find initial points for x1 and x2?? All we know right now is that they need to encompass the max?  One needs pos deriv and one needs neg
#'      method 1: gen random number and calculate h_prime
#'      method 2: find 2 stdevs from mean, check if they fit criteria
#'      
#'  3. complete s(x)
#'  4. need method to accept or reject and update ( both outputs and z )
#'  
#'          
#' @name Cadapt_reject_sample 
#' @rdname adapt_reject_sample
#' @aliases Cadapt_reject_sample
#' @exportClass Cadapt_reject_sample


library(methods) 
setClass( "Cadapt_reject_sample", 
          representation( n = "numeric", h_x = "formula", bounds = "numeric" , output = "vector", h_at_x = "vector", hprime_at_x = "vector", z = "vector", samples = "vector", x = "vector", weights = "vector", normalized_factor = "numeric" ), 
          prototype=prototype( n=50L, h_x = formula(y~-1/(sqrt(2*pi))*exp(x^2)), bounds=c(-20, 20) ) 
)

# Log normal distribution is prototype

######################################
######################################

#' Cadapt_reject_sample initialization
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("initialize", "Cadapt_reject_sample", function(.Object, n , h_x, bounds) {
  .Object@n <- n
  # Input function
  .Object@h_x <- h_x
  # Bounds of function
  .Object@bounds <- bounds
  
  # Values not input by the user
  .Object@output <- vector()
  # H(x) and H'(x) evaluated at a few points
  .Object@h_at_x <- vector()
  .Object@hprime_at_x <- vector()
  # abscissa of all points
  .Object@z <- vector()
  #  Random number for adapt/reject
  .Object@samples <- vector()
  #the weight of each piece of the integration of the upper function
  .Object@weights <- vector()
  #the x points that we have evaluated h_x for
  .Object@x <- vector()
  #the integral of the upper bound function
  .Object@normalized_factor <- numeric()
  # determine x1 and x2, draw random numbers and then determine if their first derivatives are pos and neg.
  validObject(.Object)
  .Object
})


######################################
######################################

#' Validity checks for S4 \code{adapt_reject_sample} object
#' @param object An \code{adapt_reject_sample} object

validity_ars <- function(object) {
  # checking for non-integer values
  if( is.integer( object@n ) == FALSE  ) { stop( "Input number of steps is not an integer" ) }
  if(  object@n <= 0  ) { stop( "Input number of steps is not greater than zero" ) }
} 

setValidity( "Cadapt_reject_sample", validity_ars )