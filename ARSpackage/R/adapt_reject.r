#' The adapt_reject class
#'
#' This file consists of the adapt_reject_sample class, the function that calls it and the functions it uses.
#'
#' It performs an adaptive rejection sampling process as proposed
#' by Wild and Gilks in 1992.
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{n}:}{Variable of class \code{"numeric"}, n, containing the number of points to sample}
#'    \item{\code{h_x}:}{Function of class \code{"function"}, containing the log(f(x)) to sample from.}
#'    \item{\code{h_prime}:}{Function of class \code{"function"}, containing the first derivative log(f(x)) to sample from.}
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
#' @name Cadapt_reject_sample 
#' @rdname adapt_reject_sample
#' @aliases Cadapt_reject_sample
#' @exportClass Cadapt_reject_sample
#' @author J. Bladen, L. Felberg, H.W. Tsao, S. Tu
#' @references \url{http://faculty.chicagobooth.edu/hedibert.lopes/teaching/ccis2010/1992GilksWild.pdf}
#' @return S4 \code{adapt_reject_sample} object; a vector containing
#'      \item{n} points sampled from the f(x) distribution

library(methods) 
setClass( "Cadapt_reject_sample", 
  representation( n = "numeric", h_x = "function", h_prime="function" ), 
  prototype=prototype( n=50L, h_x = function(x,mu=0, sigma=1){-1/(2*sigma^2)*(x-mu)^2}, 
                       h_prime = function(x,mu=0, sigma=1){-1/sigma^2*(x-mu)} ) )

validity_ars <- function(object) {
   # checking for non-integer values
  if( is.integer( n ) == FALSE  ) { stop( "Input number of steps is not an integer" ) }
  if(  n <= 0  ) { stop( "Input number of steps is not greater than zero" ) }
}
setValidity("Cadapt_reject_sample", validity_ars)

######################################
######################################

#' Random generating first two points
#' @param object An object

setGeneric("gen_x", function(object){standardGeneric("gen_x")})


######################################
######################################

#' Cadapt_reject_sample generating first two points
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("gen_x", signature = "Cadapt_reject_sample", function(object) {
    return( 1 )
    
} )

######################################
######################################

#' Cadapt_reject_sample initialization
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("initialize", "Cadapt_reject_sample", function(.Object, n , h_x, h_prime) {
     .Object@n <- n
     .Object@output <- list()
     .Object@h_x <- h_x
     .Object@h_prime <- h_prime
     # determine x1 and x2, draw random numbers and then determine if their first derivatives are pos and neg.
     .Object@x <- .Object@gen_x( )
})
        
######################################
######################################

#' Cadapt_reject_sample show
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("show", signature = "Cadapt_reject_sample", function(object) {

} )

######################################
######################################

#' Error checking generic
#' @param object An object

setGeneric("error_check", function(object){standardGeneric("error_check")})


######################################
######################################

#' Cadapt_reject_sample error_check
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods


setMethod("error_check", signature = "Cadapt_reject_sample", function(object) {
    return( 1 )
} )

######################################
######################################

#' Upper generic
#' @param object An object

setGeneric("upper", function(object){standardGeneric("upper")})


######################################
######################################

#' Cadapt_reject_sample upper
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods


setMethod("upper", signature = "Cadapt_reject_sample", function(object) {
    return( 1 )  
} )

######################################
######################################

#' Lower generic
#' @param object An object

setGeneric("lower", function(object, x_st, ... ){standardGeneric("lower")})


######################################
######################################

#' Cadapt_reject_sample lower
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("lower", signature = "Cadapt_reject_sample", function(object, x_star) {
  # find where x_star is in the range
  m <- as.integer( x_star > object@x )
  j <- sum( n )
  j_plus_one <- j + 1
  l_x_star <- (( object@x[j_plus_one] - x_star)*object@h_x(j) + (x_star - object@x[j])*object@h_x(j_plus_one) ) / ( object@x(j_plus_one) - object@x(j) ) 
  return( l_x_star )
} )

######################################
######################################

#' Sanple generic
#' @param object An object

setGeneric("sample", function(object){standardGeneric("sample")})


######################################
######################################

#' Cadapt_reject_sample sample
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods


setMethod("sample", signature = "Cadapt_reject_sample", function(object) {
    return( 1 )  
} )

######################################
######################################

#' Update generic
#' @param object An object

setGeneric("update", function(object){standardGeneric("update")})


######################################
######################################

#' Cadapt_reject_sample update
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("update", signature = "Cadapt_reject_sample", function(object) {
    return( 1 )  
} )

######################################
######################################

#' The adapt_reject function
#'
#' This calls the class Cadapt_reject_sample and it's methods.
#'
#' @param n_samples Number of samples desired from distribution
#' @param log_fx Log of function to sample from
#' @param log_fx_prime First derivative of log of function to sample from
#' @return S4 \code{adapt_reject_sample} object; a vector containing
#'      \item{n} points sampled from the f(x) distribution

a_r_s <- function( n_samples, log_fx, log_fx_prime, ... ){
    
    # Initialize new ARS class
    ars_class <- new( Cadapt_reject_sample, n=n_samples, h_x = log_fx, h_prime = log_fx_prime )
    
    # While we do not have enough n samples, continue to sample
    while( length( .Object@output < n ) ){
        .Object@u_x <- .Object@upper()
        .Object@s_x <- .Object@s()
        .Object@samples <- .Object@sample()
        .Object@u_x_star <- .Object@upper()
        .Object@u_x_star <- .Object@lower()
        .Object@acc_rej()
        .Object@update
    }
    
    return( ars_class )
}