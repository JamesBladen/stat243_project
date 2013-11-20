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
     .Object@output <- vector()
     .Object@h_x <- h_x
     .Object@h_prime <- h_prime
     .Object@z <- vector()
     # determine x1 and x2, draw random numbers and then determine if their first derivatives are pos and neg.
     .Object@x <- .Object@gen_x( )
})
        
######################################
######################################

#' Cadapt_reject_sample show
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("show", signature = "Cadapt_reject_sample", function(object) {
    print(" The number of samples taken:")
    print(object@n)
    print(" The samples taken:")
    print(object@output)
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

#' S(x) generic
#' @param object An object

setGeneric("s_x", function(object){standardGeneric("s_x")})


######################################
######################################

#' Cadapt_reject_sample s(x)
#' 
#' Function to normalize the upper bounds of log(f(x))
#' 
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods
#' 

setMethod("s_x", signature = "Cadapt_reject_sample", function(object){
    #Normalizing u_x
    #The basic idea is just to find out the integration of u_x on the domain by adding the integration of each piecewise of u_x
    #Seprate the domain into 3 parts: [x_0,z[1]],[z[1],z[k-1]],[z[k-1],x_a],x_0 and x_a are lower and upper bounds of the domain. If domain is R, then x_0=-inf, x_a=inf
    #Use a forloop to calculate the integrations in [z[1],z[k-1]],each piece i is a line with slope h_prime(x[i]) and pass through the point (x[i],h(x[i])), and is from z[i-1] to z[i] 
    k <- length(object@x)
    for (i in 2:(k-1)){
        piecewise_function<-function( x_prime, h_x_i, h_xprime_i, x_i ){
            exp(h_xprime_i*(x_prime-x_i)+ h_x_i )
        }
    
        hx_i <- object@h_x(object@x[i])
        hx_p_i <- object@h_prime(object@x[i])
        xi <- object@x[i]
        
        integration[i-1]<-integrate(piecewise_function,z[i-1],z[i], h_x_i=hx_i, h_xprime_i=hx_p_i, x_i=xi)[[1]]
    }
    #The function for u_x in [x_0,z[1]]
    piecewise_function_1<-function(x_prime, h_x_1, h_xprime_1, x_1){
        exp(object@h_prime(object@x[1])*(x_prime-object@x[1])+object@h_x(object@x[1]))
    }
    
    #The function for u_x in [z[k],x_a]
    piecewise_function_k<-function(x_prime, h_x_k, h_xprime_k, x_k){
        exp(object@h_prime(object@x[k])*(x_prime-object@x[k])+object@h_x(object@x[k]))
    }
    
    # Calculating constants for integration    
    hx_1 <- object@h_prime(object@x[1])
    hxprime_1 <- object@h_x(object@x[1])
    x1 <- object@x[1]
    hx_k <- object@h_x(object@x[k])
    hxprime_k <- object@h_prime(object@x[k])
    xk <- object@x[k]
    
    #Add up all 3 parts
    normalized_factor<-integrate(piecewise_function_1,x_0,z[1], hx_1, hxprime_1, x1)+integrate(piecewise_function_k,z[k-1],x_a, hx_k, hxprime_k, xk)+sum(integration)
    
    weights <- 0
})

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
    #Figure out the vector of z
    z<-vector()
    # k is number of xs
    for (i in 1:(k-1)){
        z[i]<-(object@h_x(object@x[i+1])-object@h_x(object@x[i])-object@x[i+1]*object@h_prime(object@x[i+1])+object@x[i]*obejct@h_prime(object@x[i]))/(object@h_prime(object@x[i])-object@h_prime(object@x[i+1]))
    }
    #Calculate u of x star using the same method as we calculate l of x star
    M<-as.integer(x_star>z)
    J<-sum(M)
    J_plus_one<-J+1
    u_x_star<-object@h_x(object@x[J_plus_one])+(x_star-object@x[J_plus_one])*object@h_prime(object@x[J_plus_one])
    return(u_x_star)
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
    samples <- vector()
    # Sample uniform random number
    samples[1] <- runif( 1, min = 0, max = 1 )
    
    # Sample x_star from sk(x)
    samples[2] <- object@sample_from_S(  )
    
    
    return( 1 )  
} )

######################################
######################################

#' Sample from S(x) generic
#' @param object An object

setGeneric("sample_from_S", function(object){standardGeneric("sample_from_S")})


######################################
######################################

#' Cadapt_reject_sample sample_from_S
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("sample_from_S", signature = "Cadapt_reject_sample", function(object) {
    
    
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