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
     .Object@weights<-vector()
     .Object@normalized_factor<-numeric()
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
    #Figure out the vector of z
    z<-vector()
    # k is number of xs
    for (i in 1:(k-1)){
        z[i+1]<-(object@h_x(object@x[i+1])-object@h_x(object@x[i])-object@x[i+1]*object@h_prime(object@x[i+1])+object@x[i]*obejct@h_prime(object@x[i]))/(object@h_prime(object@x[i])-object@h_prime(object@x[i+1]))
    }
    z[1]<-x_0
    z[k+1]<-x_a
    for (i in 1:k){
        piecewise_function<-function( x_prime, h_x_i, h_xprime_i, x_i ){
            exp(h_xprime_i*(x_prime-x_i)+ h_x_i )
        }
    
        hx_i <- object@h_x(object@x[i])
        hx_p_i <- object@h_prime(object@x[i])
        xi <- object@x[i]
        
        integration[i]<-integrate(piecewise_function,z[i],z[i+1], h_x_i=hx_i, h_xprime_i=hx_p_i, x_i=xi)[[1]]
    }   
        normalized_factor<-sum(integration) 
    
        weights<-integration/normalized_factor
    
    # Calculating constants for integration    
    hx_1 <- object@h_prime(object@x[1])
    hxprime_1 <- object@h_x(object@x[1])
    x1 <- object@x[1]
    hx_k <- object@h_x(object@x[k])
    hxprime_k <- object@h_prime(object@x[k])
    xk <- object@x[k]
    
    
    object@weights<-weights
    object@normalized_factor<-normalized_factor
    object@z<-z
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
    k<-length(object@x)
    region_x_star<-sample(1:k,1,prob=object@weights)
    a<-object@h_prime(object@x[region_x_star])
    b<-object@h(object@x[region_x_star])-object@h_prime(object@x[region_x_star])*object@x[region_x_star]
    inverse_CDF<-function(x_prime){
        (log(a*x_prime/(object@normalized_factor*exp(b))+exp(a*object@z[region_x_star])))/a
    }
    sample_uniform<-runif(1)
    x_star<-inverse_CDF(sample_uniform) 
    samples[2] <- x_star
    
    
    return(samples)  
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

