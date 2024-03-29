#####################################################
#####################################################

test <- function( nsamp, function_x, bds, guess_mode=-999 ){
    
    
    #####################################################
    #####################################################
    # We implemented the following tests for each method.  There are a few different variants of the code.  It can be run on different densities and/or it can be run with a variety of bounds
    
    # First, we run through each method,
    
    ars_class <- new( "Cadapt_reject_sample", n=nsamp, f_x = function_x, bounds=bds )
    
    # Testing gen_x # Implentation with different bounds, should result in proper generation of starting values.  
    
    # There are four cases that can be selected no bounds, a left bound only, a right bound only and both sides bounded.  For example, the First case, no bounds is given as follows:

    # It can be tested with the following lines:
    ars_class<-gen_x(ars_class)
    ars_class@x
    ars_class@h_at_x
    ars_class@hprime_at_x
    
    # We expect to see that reasonable starting values to be chosen.  Two points with different derivatives.  X may be larger than two values, but we need atleast two points with different derivative signs.
    
    #####################################################
    #####################################################
    
    # the next method to test is the s_x function, which is computed using an inverse CDF.  This can be called using the following:
    
    ars_class<-s_x(ars_class)
    
    # Z should contain a vector of values length( ars_class@x - 1 ) + 2.  It contains the abcissa of the x vector and the user input bounds.  The values should all be between the points in x, with exception of the bounds.  
    ars_class@z
    
    # The integration of each interval of s_x, all values should be positive.
    ars_class@piecewise_integration
    
    # The total integration of s_x
    ars_class@normalized_factor
    sum(ars_class@piecewise_integration) == ars_class@normalized_factor # should be TRUE
    
    # This slot contains a vector of weights of each piece of the piecewise integration function
    ars_class@weights
    sum(ars_class@weights) == 1
    
    
    #####################################################
    #####################################################
    
    # The next function tested is sampling.  It generates a sample from U(x)
    ars_class<-sampling(ars_class)
    
    # It should return a single value
    ars_class@samples
    
    
    #####################################################
    #####################################################
    
    # Update progresses the ARS method.  It calls upper and lower and performs the ratio test.  During an update call, it compares a sampled x values to the ratio given in the paper.  It will either add that accepted number to the output, or if rejected it will compare it to the upper bounds and if accepted add it to x as well as the output, otherwise reject.    
    ars_class<-update(ars_class)
    
    # Depending on where the sample falls in the ACC/REJ strata, it should 1. just increment to @output, 2. increment all variable below or 3. reject and change nothing.
    ars_class@x
    ars_class@h_at_x
    ars_class@hprime_at_x
    ars_class@output
    
}