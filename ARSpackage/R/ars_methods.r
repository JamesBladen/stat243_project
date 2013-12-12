######################################
######################################

setGeneric("gen_x", function(object){standardGeneric("gen_x")})

#' gen_x
#' 
#' Cadapt_reject_sample method for generating first two points.  If the distribution is unbounded, then find the function's mode and pick points surrounding it.  If it's bounded on one side, we use the bound given and search until we find a point that corresponds to the opposite end of the domain with respect to their derivatives.  If bounded on both sides, use given bounds.  
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars_method



setMethod("gen_x", signature = "Cadapt_reject_sample", function(object) {
    
    
    # First case, the function is bounded on neither side 
    if(object@bounds[1]==-Inf && object@bounds[2]==Inf){
        
        # If user specified a mode, we calculate the true mode and start from (true mode-0.1). If not we start from -0.1
        if(object@guess_of_mode!=-999){
            bound<-c(object@guess_of_mode-100, object@guess_of_mode+100)
            the_mode<-optimize(object@f_x, interval=bound, maximum=T)$maximum
        }else{
            the_mode<-0
        }
        
        
        object@x[1]<-the_mode - .1
        
        h_x<-log(object@f_x(object@x[1]))
        deriv1 <- (1/object@f_x(object@x[1])) *genD(object@f_x,object@x[1])$D[1]
        object@h_at_x[1]<-h_x
        object@hprime_at_x[1]<-deriv1
        
        # If the initial point has positive slope, we move left by 1/2 each time until we get a point with negative slope and use that point as our second initial point
        if(sign(deriv1==1)){
            i<-2
            while(sign(deriv1)!=-1){
                object@x[i]<- object@x[1]+1/2*(i-1)
                h_x<-log(object@f_x(object@x[i]))
                deriv1 <- (1/object@f_x(object@x[i])) *genD(object@f_x,object@x[i])$D[1]
                
                object@h_at_x[i]<-h_x
                object@hprime_at_x[i]<-deriv1
                i<-i+1
            }
        # If the initial point has negative slope, we move right by 1/2 each time until we get a point with positive slope and use that point as our second initial point
        }else if(sign(deriv1)==-1){
            i<-2
            while(sign(deriv1)!=1){
                object@x[i]<- object@x[1]-1/2*(i-1)
                h_x<-log(object@f_x(object@x[i]))
                deriv1 <- (1/object@f_x(object@x[i])) *genD(object@f_x,object@x[i])$D[1]
                
                object@h_at_x[i]<-h_x
                object@hprime_at_x[i]<-deriv1
                i<-i+1
            }
        # If the initial point has 0 slope, then we pick one point 1/2 on the left and one point 1/2 on the right
        }else{
            object@x[2]<- object@x[1]-1/2
            h_x<-log(object@f_x(object@x[2]))
            deriv1 <- (1/object@f_x(object@x[2])) *genD(object@f_x,object@x[2])$D[1]
            object@h_at_x[2]<-h_x
            object@hprime_at_x[2]<-deriv1
            
            object@x[3]<- object@x[1]+1/2
            h_x<-log(object@f_x(object@x[3]))
            deriv1 <- (1/object@f_x(object@x[3])) *genD(object@f_x,object@x[3])$D[1]
            object@h_at_x[3]<-h_x
            object@hprime_at_x[3]<-deriv1
        }
    # Second case, the function is bounded only on the right
    }else if (object@bounds[1]==-Inf && object@bounds[2] !=Inf){
        
        # If user specified a mode, we calculate the true mode and start from (true mode+0.99), or (right bound-0.01) if (true mode+0.99) is outside the right bound. If not we start from (right bound-0.01)
        if(object@guess_of_mode!=-999){
            bound<-c(object@guess_of_mode-100, object@guess_of_mode+100)
            the_mode<-optimize(object@f_x, interval=bound, maximum=T)$maximum
            
            if(the_mode > object@bounds[2])
            {
                the_mode<-object@bounds[2]
            }
        }else{
            the_mode<-object@bounds[2]
        }
        
        
        
        object@x[1]<-the_mode + min(0.1, abs(the_mode - object@bounds[1])) -.01
        
        
        
        
        
        h_x<-log(object@f_x(object@x[1]))
        deriv1 <- (1/object@f_x(object@x[1])) *genD(object@f_x,object@x[1])$D[1]
        
        object@h_at_x[1]<-h_x
        object@hprime_at_x[1]<-deriv1
        # If the initial point has negative slope, we move left by 1/2 each time until we get a point with positive slope and use that point for second initial point 
        if(sign(deriv1)==-1){
            i<-2
            while(sign(deriv1)!=1){
                object@x[i]<- object@x[1]- 1/2*(i-1)
                h_x<-log(object@f_x(object@x[i]))
                deriv1 <- (1/object@f_x(object@x[i])) *genD(object@f_x,object@x[i])$D[1]
                
                object@h_at_x[i]<-h_x
                object@hprime_at_x[i]<-deriv1
                i<-i+1
            }
        # If the initial point has positive slope this means the function is monotonically increasing, then we take a point which is (initial point-1/2) for our second initial point
        }else{
            object@x[2]<-object@x[1] -1/2
            h_x<-log(object@f_x(object@x[2]))
            deriv1 <- (1/object@f_x(object@x[2])) *genD(object@f_x,object@x[2])$D[1]
            object@h_at_x[2]<-h_x
            object@hprime_at_x[2]<-deriv1
        }
    # Third case, the function is bounded only on the left
    }else if (object@bounds[1]!=-Inf && object@bounds[2] ==Inf){
        
        
        
      # If user specified a mode, we calculate the true mode and start from (true mode-0.99), or (left bound+0.01) if (true mode-0.99) is outside the left bound. If not we start from (left bound+0.01)
        if(object@guess_of_mode!=-999){
            bound<-c(object@guess_of_mode-100, object@guess_of_mode+100)
            the_mode<-optimize(object@f_x, interval=bound, maximum=T)$maximum
            
            if(the_mode < object@bounds[1])
            {
                the_mode<-object@bounds[1]
            }
        }else{
            the_mode<-object@bounds[1]
        }
        
        
        
        
        
        object@x[1]<-the_mode - min(0.1, abs(the_mode - object@bounds[1])) +.01
        
        
        
        
        h_x<-log(object@f_x(object@x[1]))
        deriv1 <- (1/object@f_x(object@x[1])) *genD(object@f_x,object@x[1])$D[1]
        
        object@h_at_x[1]<-h_x
        object@hprime_at_x[1]<-deriv1
        # If the initial point has positive slope, we move right by 1/2 each time until we get a point with negative slope and use that point for second initial point 
        if(sign(deriv1)==1)
        {
            i<-2
            while(sign(deriv1)!=-1){
                object@x[i]<- object@x[1]+ 1/2*(i-1)
                h_x<-log(object@f_x(object@x[i]))
                deriv1 <- (1/object@f_x(object@x[i])) *genD(object@f_x,object@x[i])$D[1]
                
                object@h_at_x[i]<-h_x
                object@hprime_at_x[i]<-deriv1
                i<-i+1
            }
        # If the initial point has negative slope this means the function is monotonically decreasing, then we take a point which is (initial point+1/2) for our second initial point
        }else{
            object@x[2]<-object@x[1] +1/2
            h_x<-log(object@f_x(object@x[2]))
            deriv1 <- (1/object@f_x(object@x[2])) *genD(object@f_x,object@x[2])$D[1]
            object@h_at_x[2]<-h_x
            object@hprime_at_x[2]<-deriv1
        }
    # Fourth case, the function is bounded on both sides. We pick one initial point to be the (left bound+0.1) and one initial point to be the (right bound-0.1)
    }else if(object@bounds[1]!=-Inf && object@bounds[2] !=Inf){
        
        object@x[1]<-object@bounds[1]+0.1
        object@x[2]<-object@bounds[2]-0.1
        h_x<-log(object@f_x(object@x[1]))
        deriv1 <- (1/object@f_x(object@x[1])) *genD(object@f_x,object@x[1])$D[1]
        object@h_at_x[1]<-h_x
        object@hprime_at_x[1]<-deriv1
        
        h_x<-log(object@f_x(object@x[2]))
        deriv1 <- (1/object@f_x(object@x[2])) *genD(object@f_x,object@x[2])$D[1]
        object@h_at_x[2]<-h_x
        object@hprime_at_x[2]<-deriv1
    }
    
    return(object)
} )



######################################
######################################

setGeneric("ev_h", function(object){standardGeneric("ev_h")})

#' ev_h
#' 
#' Cadapt_reject_sample method that evaluates the log(f(x)) for a given x and the derivative as well. 
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars_method


setMethod("ev_h", signature = "Cadapt_reject_sample", function(object) {
    
    
    h_x<-log(object@f_x(object@samples[2]))
    h_value <- (1/object@f_x(object@samples[2])) *genD(object@f_x,object@samples[2])$D[1]
    
    
    return( c(h_x,h_value) )
} )


######################################
######################################

setGeneric("s_x", function(object){standardGeneric("s_x")})

#' s(x)
#' 
#' Cadapt_reject_sample method to normalize the upper bounds of log(f(x)).  Multiple objective are performed here. The most important being the calculation of the abcissa vector Z.  Additionally, the weights and exact values of the piecewise integration of each interval and the normalization factor for the entire upper bound are calculated and the x's, their evaluations and their derivatives are sorted by x.  
#' 
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars_method


setMethod("s_x", signature = "Cadapt_reject_sample", function(object){
    k <- length(object@x)
    z<-vector()
    # Sort x, h_at_x and hprime_at_x by x
    mat<-cbind(object@x,object@h_at_x,object@hprime_at_x)
    mat_sort<-mat[order(mat[,1]),]
    # Create lower and upper bound vector for sorted x, h_at_x and hprime_at_x
    x_low<-mat_sort[,1][-k]
    x_high<-mat_sort[,1][-1]
    h_at_x_low<-mat_sort[,2][-k]
    h_at_x_high<-mat_sort[,2][-1]
    hprime_at_x_low<-mat_sort[,3][-k]
    hprime_at_x_high<-mat_sort[,3][-1]
    # Create abcissa vector z
    z_prime<-(h_at_x_high-h_at_x_low-x_high*hprime_at_x_high+x_low*hprime_at_x_low)/(hprime_at_x_low-hprime_at_x_high)
    z<-c(object@bounds[1],z_prime,object@bounds[2])
    # Create lower and upper bound vector of z
    z_low<-z[-(k+1)]
    z_high<-z[-1]
    # Calculate the integration of exp(u_x) in each interval. If hprime_at_x is 0 so the integration is NaN, we calculate it manually 
    piecewise_integration<-(1/mat_sort[,3])*(exp(mat_sort[,3]*z_high+mat_sort[,2]-mat_sort[,3]*mat_sort[,1])-(exp(mat_sort[,3]*z_low+mat_sort[,2]-mat_sort[,3]*mat_sort[,1])))
    piecewise_integration[which(piecewise_integration=="NaN")]<-exp(mat_sort[,2][which(piecewise_integration=="NaN")]) *(z_high[which(piecewise_integration=="NaN")]-z_low[which(piecewise_integration=="NaN")])
    # Sum up the integration on each interval to find the normalizing factor for exp(u_x), and find the weights for integration on each interval
    normalized_factor<-sum(piecewise_integration)
    weights<-piecewise_integration/normalized_factor
    
    #check for log cancavity
    
    if(sum(weights<0)!=0)
    {
        print("function is not log concave")
        break
    }
    
    object@weights<-weights
    object@normalized_factor<-normalized_factor
    object@z<-z
    object@mat_sorted<-mat_sort
    object@piecewise_integration<-piecewise_integration
    return(object)
})


######################################
######################################

setGeneric("sampling", function(object){standardGeneric("sampling")})

#' sampling
#' 
#' Method to sample from s_x.  The basic algorithm is as follows: 1. Determine an interval to sample from using the weights of integration of the function on each interval, computed in the s_x method. 2.  Use inverse CDF method to sample from within the selected interval.  Return the object with new sample.
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars_method


setMethod("sampling", signature = "Cadapt_reject_sample", function(object) {
    samples <- vector()
    # Sample uniform random number
    object@samples[1] <- runif( 1, min = 0, max = 1 )
    
    # Sample x_star from s_k(x) 
    k<-length(object@x)
    # First sample the interval in which we will sample x_star
    region_x_star<-sample(1:k,1,prob=object@weights)
    # Using inverse CDF to sample x_star within the interval. Seperate cases for hprime_at_x is 0 and is not 0
    a<-object@mat_sorted[,3][region_x_star]
    
    if(a==0){
        x_star<-runif(1,object@z[region_x_star],object@z[region_x_star +1])
    }else{
        b<-object@mat_sorted[,2][region_x_star]-object@mat_sorted[,3][region_x_star]*object@mat_sorted[,1][region_x_star]
        inverse_CDF<-function(x_prime){
            (log(a*x_prime*object@piecewise_integration[region_x_star]/exp(b)+exp(a*object@z[region_x_star])))/a
        }
        sample_uniform<-runif(1)
        x_star<-inverse_CDF(sample_uniform) 
    }
    
    object@samples[2] <- x_star
    
    return( object )  
} )


######################################
######################################

setGeneric("upper", function(object){standardGeneric("upper")})

#' upper
#' 
#' Cadapt_reject_sample method to evaluate the upper bound of x_star.
#' 
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars_method



setMethod("upper", signature = "Cadapt_reject_sample", function(object) {
    x_star<-object@samples[2]
    # Find which interval of z does x_star fall in 
    M<-as.integer(x_star > object@z)
    J<-sum(M)
    # Evaluate u_j(x_star)
    u_x_star<-object@mat_sorted[,2][J]+(x_star-object@mat_sorted[,1][J])*object@mat_sorted[,3][J]
    return(u_x_star)
} )


######################################
######################################

setGeneric("lower", function(object, x_st, ... ){standardGeneric("lower")})

#' lower
#' 
#' Cadapt_reject_sample method to evaluate the lower bound of x_star.
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars_method


setMethod("lower", signature = "Cadapt_reject_sample", function(object) {
    x_star<-object@samples[2]
    # Find which interval of x does x_star fall in 
    m <- as.integer( x_star > object@x )
    j <- sum( m )
    j_plus_one <- j + 1
    # Evaluate l_j(x_star)
    l_x_star <- (( object@mat_sorted[,1][j_plus_one] - x_star)*object@mat_sorted[,2][j] + (x_star- object@mat_sorted[,1][j])*object@mat_sorted[,2][j_plus_one] ) / ( object@mat_sorted[,1][j_plus_one] - object@mat_sorted[,1][j] ) 
    return( l_x_star )
} )


######################################
######################################

setGeneric("update", function(object){standardGeneric("update")})

#' update
#'
#' Cadapt_reject_sample method to determine which ACC/REJ criteria a given sampled value fits into and updates the samples and x values accordingly.  
#' 
#' 
#' @param object \code{\linkS4class{Cadapt_reject_sample}} object
#' @rdname ars_method



setMethod("update", signature = "Cadapt_reject_sample", function(object) {
    w<-object@samples[1]
    
    if(sum(object@samples[2]< object@x)==0  || sum(object@samples[2]< object@x)==length(object@x))
    {
        ratio<--1
    }else{
        ratio<-exp(lower(object) - upper(object))
    }
    
    #calculate the ratio of lower to upper  
    if(w<=ratio){
        #if we are within this first ratio, add to output
        object@output<-c(object@output,object@samples[2])
    }else{
        #if we aren't in the first ratio, calc hstar and hprimestar
        hvals <- ev_h(object)
        hstar <- hvals[1]
        hprimestar <- hvals[2]
        
        
        
        
        ratio<-exp(hstar-upper(object))
        
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
    
    return(object)
} )
