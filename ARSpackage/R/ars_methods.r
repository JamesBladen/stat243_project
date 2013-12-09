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

#' Random generating first two points
#' @param object An object

setGeneric("gen_x", function(object){standardGeneric("gen_x")})


######################################
######################################

#' Cadapt_reject_sample generating first two points
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("gen_x", signature = "Cadapt_reject_sample", function(object) {
  object@x <- object@bounds
  
  
  
  
  h_x<-log(object@f_x(object@x[1]))
  h_value <- (1/object@f_x(object@x[1])) *genD(object@f_x,object@x[1])$D[1]
  
  
  object@h_at_x[1]<-h_x
  object@hprime_at_x[1]<-h_value
  
  h_x<-log(object@f_x(object@x[2]))
  h_value <- (1/object@f_x(object@x[2])) *genD(object@f_x,object@x[2])$D[1]
  
  object@h_at_x[2]<-h_x
  object@hprime_at_x[2]<-h_value
  return(object)
} )



######################################
######################################

#' eval_h checking generic
#' @param object An object

setGeneric("ev_h", function(object){standardGeneric("ev_h")})


######################################
######################################

#' Cadapt_reject_sample eval_h
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("ev_h", signature = "Cadapt_reject_sample", function(object) {
  
  
  h_x<-log(object@f_x(object@samples[2]))
  h_value <- (1/object@f_x(object@samples[2])) *genD(object@f_x,object@samples[2])$D[1]
  
  
  return( c(h_x,h_value) )
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
  # k is number of xs 
  k <- length(object@x)
  #Figure out the vector of z
  z<-vector()
  
  mat<-cbind(object@x,object@h_at_x,object@hprime_at_x)
  mat_sort<-mat[order(mat[,1]),]
  x_low<-mat_sort[,1][-k]
  x_high<-mat_sort[,1][-1]
  h_at_x_low<-mat_sort[,2][-k]
  h_at_x_high<-mat_sort[,2][-1]
  hprime_at_x_low<-mat_sort[,3][-k]
  hprime_at_x_high<-mat_sort[,3][-1]
  z_prime<-(h_at_x_high-h_at_x_low-x_high*hprime_at_x_high+x_low*hprime_at_x_low)/(hprime_at_x_low-hprime_at_x_high)
  z<-c(object@bounds[1],z_prime,object@bounds[2])
  
  z_low<-z[-(k+1)]
  z_high<-z[-1]
  piecewise_integration<-(1/mat_sort[,3])*(exp(mat_sort[,3]*z_high+mat_sort[,2]-mat_sort[,3]*mat_sort[,1])-(exp(mat_sort[,3]*z_low+mat_sort[,2]-mat_sort[,3]*mat_sort[,1])))
  normalized_factor<-sum(piecewise_integration)
  weights<-piecewise_integration/normalized_factor
  
  object@weights<-weights
  object@normalized_factor<-normalized_factor
  object@z<-z
  object@mat_sorted<-mat_sort
  return(object)
})


######################################
######################################

#' Sample generic
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
  object@samples[1] <- runif( 1, min = 0, max = 1 )
  
  # Sample x_star from sk(x)
  k<-length(object@x)
  region_x_star<-sample(1:k,1,prob=object@weights)
  a<-object@hprime_at_x[region_x_star]
  b<-object@h_at_x[region_x_star]-object@hprime_at_x[region_x_star]*object@x[region_x_star]
  inverse_CDF<-function(x_prime){
    (log(a*x_prime/(object@normalized_factor*exp(b))+exp(a*object@z[region_x_star])))/a
  }
  sample_uniform<-runif(1)
  x_star<-inverse_CDF(sample_uniform) 
  object@samples[2] <- x_star
  
  return( object )  
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

#' Update generic
#' @param object An object

setGeneric("update", function(object){standardGeneric("update")})


######################################
######################################

#' Cadapt_reject_sample update
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars-methods

setMethod("update", signature = "Cadapt_reject_sample", function(object) {
  w<-object@samples[1]
  
  #calculate the ratio of lower to upper
  ratio<-exp(object@lower( object, object@samples[2] ) - object@upper( object, object@samples[2]  ) )
  
  if(w<=ratio){
    #if we are within this first ratio, add to output
    object@output<-c(object@output,object@samples[2])
  }else{
    #if we aren't in the first ratio, calc hstar and hprimestar
    hvals <- object@eval_h()
    hstar <- hvals[1]
    hprimestar <- hvals[2]
    
    ratio<-exp( hstar-object@upper( object, object@samples[2]  ) )
    
    if(w<=ratio){
      #if we accept, update our x, put xstar in output, and add hstar/hprimestar 
      object@output<-c(object@output,object@samples[2])
      object@x<-c(object@x,object@samples[2])
      object@h_at_x<-c(object@h_at_x,hstar)
      object@hprime_at_x<-c(object@hprime_at_x,hprimestar)
    }else{
      #update our x and add hstar hprimexstar
      object@x<-c(object@x,object@samples[2])
      object@h_at_x<-c(object@h_at_x,hstar)
      object@hprime_at_x<-c(object@hprime_at_x,hprimestar)
    }
  } 
} )
