#' The EXAMPLE class
#'
#' This class contains an example. This line goes into the description
#'
#' This line and the next ones go into the details.
#' This line thus appears in the details as well.
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{Matrix of class \code{"numeric"}, containing data from slot1}
#'    \item{\code{slot2}:}{Object of class \code{"character"}, containing data that needs to go in slot2.}
#'  }
#'
#' @note You can still add notes
#' @name EXAMPLE 
#' @rdname EXAMPLE
#' @aliases EXAMPLE-class
#' @exportClass EXAMPLE
#' @author Joris Meys

# 1. Initialize
#   i) x1, x2
#   ii) inputs: h(x) and h'(x), n (number of points to sample), optional: domain etc
#   iii) error checks: make sure that the function is concave up and the function lies within U(x) and L(x).  
#         Check that x1 has a positive slope and X2 has a negative slope. Check that the sample size is positive and an integer.
# 2) Objects/methods:
#   i) U(x) and S(x): z(x), equations for tangent lines
#   ii) List of x points
#   iii) list of sampled points
#   iv) l(x)
#   v) sample function from s(x) and uniform random number
#   vi) update steps
#   vii) error checking

library(methods) 
setClass( "adapt_reject_sample", 
  representation( n = "numeric", h_x = "function", h_prime="function" ), 
  prototype=prototype( n_steps=50L, h_x = function(x,mu=0, sigma=1){-1/(2*sigma^2)*(x-mu)^2}, 
                       h_prime = function(x,mu=0, sigma=1){-1/sigma^2*(x-mu)} ) )

validity_ars <- function(object) {
   # checking for non-integer values
  if( is.integer( n ) == FALSE  ) { stop( "Input number of steps is not an integer" ) }
  if(  n <= 0  ) { stop( "Input number of steps is not greater than zero" ) }
}
setValidity("adapt_reject_sample", validity_ars)

setMethod("initialize", "adapt_reject_sample", function(.Object, n , h_x, h_prime) {
     .Object@n <- n
     .Object@output <- list()
     .Object@h_x <- h_x
     .Object@h_prime <- h_prime
     # determine x1 and x2, draw random numbers and then determine if their first derivatives are pos and neg.
     .Object@x <- .Object@gen_x( )
     while( length( .Object@output < n ) ){
       .Object@u_x <- .Object@upper()
       .Object@s_x <- .Object@s()
       .Object@samples <- .Object@sample()
       .Object@u_x_star <- .Object@upper()
       .Object@u_x_star <- .Object@lower()
       .Object@acc_rej()
       .Object@update
     }
       
})
        
                      

setMethod("show", signature = "adapt_reject_sample", function(object) {

} )

setMethod("gen_x", signature = "adapt_reject_sample", function(object) {
  
} )

setMethod("error_check", signature = "adapt_reject_sample", function(object) {
  
} )

setMethod("upper", signature = "adapt_reject_sample", function(object) {
  
} )

setMethod("lower", signature = "adapt_reject_sample", function(object, x_star) {
  # find where x_star is in the range
  m <- as.integer( x_star > object@x )
  j <- sum( n )
  j_plus_one <- j + 1
  l_x_star <- (( object@x[j_plus_one] - x_star)*object@h_x(j) + (x_star - object@x[j])*object@h_x(j_plus_one) ) 
                / ( object@x(j_plus_one) - object@x(j) ) 
  return( l_x_star )
} )

setMethod("sample", signature = "adapt_reject_sample", function(object) {
  
} )

setMethod("update", signature = "adapt_reject_sample", function(object) {
  
} )