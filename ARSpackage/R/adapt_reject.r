######################################
######################################

#' The adapt_reject function
#'
#' This calls the class Cadapt_reject_sample and its methods.  The vector of samples is accessible via \var{ans output}.
#'
#' @param n_samples: Number of samples desired from distribution
#' @param fx: Function to sample from
#' @param bounds: Bounds of function of interest.  The default is an unbounded function
#' @return S4 \code{adapt_reject_sample} object; a vector containing \code{n} points sampled from the f(x) distribution
#' 

a_r_s <- function( n_samples, fx, bounds=c(-Inf, Inf), ... ){
  
  # Initialize new ARS class
  ars_class <- new( "Cadapt_reject_sample", n=n_samples, f_x = fx, bounds, ... )
  ars_class <- gen_x( ars_class )
  
  #print( ars_class@output )
  
  # While we do not have enough n samples, continue to sample
  while( length( ars_class@output) < n_samples ){
    
    ars_class <- s_x( ars_class )
    ars_class <- sampling( ars_class )
    ars_class <- update( ars_class )
    
  }

  return(  ars_class@output )
}

######################################
######################################
#' The adapt_reject class
#'
#' This class contains all the methods used to perform an AR sampling.  
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{n}:}{Variable of class \code{"numeric"}, n, containing the number of points to sample, taken as user input.}
#'    \item{\code{f_x}:}{Function of class \code{"function"}, containing the f(x) to sample from, taken as user input.}
#'    \item{\code{bounds}:}{Variable of class \code{"numeric"}, n, containing the bounds of the function, taken as user input.}
#'    \item{\code{output}:}{Variable of class \code{"vector"}, containing sampled points to return to user.}
#'    \item{\code{h_at_x}:}{Variable of class \code{"vector"}, containing computed log(f(x)) values at all x values}
#'    \item{\code{hprime_at_x}:}{Variable of class \code{"vector"}, containing computed derivative of log(f(x)) values at all x values}

#'    \item{\code{z}:}{Variable of class \code{"vector"}, containing abscissae of upper bound function.}
#'    \item{\code{samples}:}{Variable of class \code{"vector"}, containing random numbers generated by s(x) and unif.}
#'    \item{\code{x}:}{Variable of class \code{"vector"}, containing x values used in ARS.}
#'    \item{\code{weights}:}{Variable of class \code{"vector"}, containing sampled points to return to user.}    
#'    \item{\code{output}:}{Variable of class \code{"numeric"}, containing sampled points to return to user.}
#'    \item{\code{mat_sorted}:}{Variable of class \code{"matrix"}, containing x values, their corresponding h and h prime values, sorted by increasing x.}
#'  }
#'          
#' @name Cadapt_reject_sample 
#' @aliases Cadapt_reject_sample
#' @exportClass Cadapt_reject_sample

library(methods) 
library(numDeriv)
setClass( "Cadapt_reject_sample", 
          representation( n = "numeric", f_x = "function", bounds = "numeric" , output = "vector", h_at_x = "vector", hprime_at_x = "vector", z = "vector", samples = "vector", x = "vector", weights = "vector", normalized_factor = "numeric", mat_sorted="matrix",piecewise_integration="vector" ), 
          prototype=prototype( n=50L, f_x = function(x){(-1/(2*1^2)*exp((x-0)^2))}, bounds=c(-20, 20) ) 
)

# Log normal distribution is prototype

######################################
######################################

#' Cadapt_reject_sample initialization: method to intialize the ARS class for sampling.  Will store values input from user and will also initialize empty arrays for all other slots.
#' 
#' 
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @param n \code{numeric} determining the number of samples to obtain
#' @param f_x \code{function} for distribution to sample from
#' @param bounds \code{vector} of distribution bounds


setMethod("initialize", "Cadapt_reject_sample", function(.Object, n , f_x, bounds) {
    # User inputs
    
    # number of samples to take
    .Object@n <- n
    # Input function
    .Object@f_x <- f_x
    # Bounds of function
    .Object@bounds <- bounds
  
  # Values not input by the user
    
    #the x points that we have evaluated h_x for
    .Object@x <- vector()
    # samples to return to user
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
  #the integral of the upper bound function
  .Object@normalized_factor <- numeric()
  #The matrix of x, h_at_x and hprime_at_x sorted based on x
  .Object@mat_sorted<-matrix()
  #The vector of the integration of u_x on each piece
  .Object@piecewise_integration<-vector()
  
  validObject(.Object)
  .Object
})


######################################
######################################

#' Validity checks for S4 \code{adapt_reject_sample} object: want to ensure at creation that the number of samples desired is a positive integer
#' @param object An \code{adapt_reject_sample} object

validity_ars <- function(object) {
  # checking for non-integer values
  if( object@n %% 1 != 0  ) { stop( "Input number of steps is not an integer" ) }
  if(  object@n <= 0  ) { stop( "Input number of steps is not greater than zero" ) }
} 

setValidity( "Cadapt_reject_sample", validity_ars )
