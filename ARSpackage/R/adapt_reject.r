
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

a_r_s <- function( n_samples, log_fx, log_fx_prime, ... ){
    
    # Initialize new ARS class
    ars_class <- new( Cadapt_reject_sample, n=n_samples, h_x = log_fx, h_prime = log_fx_prime )
    
    # While we do not have enough n samples, continue to sample
    while( length( ars_class@output < n ) ){
        ars_class@u_x <- ars_class@upper()
        ars_class@s_x <- ars_class@s()
        ars_class@samples <- ars_class@sample()
        ars_class@u_x_star <- ars_class@upper()
        ars_class@u_x_star <- ars_class@lower()
        ars_class@acc_rej()
        ars_class@update
    }
    
    return( ars_class )
}


#' The adapt_reject class
#'
#' This class contains all the methods used to perform an AR sampling.  
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{n}:}{Variable of class \code{"numeric"}, n, containing the number of points to sample}
#'    \item{\code{h_x}:}{Function of class \code{"function"}, containing the log(f(x)) to sample from.}
#'    \item{\code{h_prime}:}{Function of class \code{"function"}, containing the first derivative log(f(x)) to sample from.}
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
  representation( n = "numeric", h_x = "function", h_prime="function" ), 
  prototype=prototype( n=50L, h_x = function(x,mu=0, sigma=1){-1/(2*sigma^2)*(x-mu)^2}, 
                       h_prime = function(x,mu=0, sigma=1){-1/sigma^2*(x-mu)} ) 
)


######################################
######################################

#' Cadapt_reject_sample initialization
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("initialize", "Cadapt_reject_sample", function(.Object, n , h_x, h_prime) {
     .Object@n <- n
     .Object@output <- vector()
     .Object@h_x <- h_x
     .Object@h_prime <- h_prime
     .Object@z <- vector()
     # determine x1 and x2, draw random numbers and then determine if their first derivatives are pos and neg.
     .Object@x <- .Object@gen_x( )
     validObject(.Object)
	.Object
})
   

######################################
######################################

#' Validity checks for S4 \code{adapt_reject_sample} object
#' @param object An \code{adapt_reject_sample} object

validity_ars <- function(object) {
   # checking for non-integer values
  if( is.integer( n ) == FALSE  ) { stop( "Input number of steps is not an integer" ) }
  if(  n <= 0  ) { stop( "Input number of steps is not greater than zero" ) }
} 

setValidity( "Cadapt_reject_sample", validity_ars )