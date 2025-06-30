
#kde approximation for the ct values 
# B number of sims 

#kde function approximation 

kde_approx<-function(B,x,viral_loads,intercept, slope ){
  ct_threshold<-40 
#  den_estimates <- NULL
  all_kdes <- list(length=51)
  
  all_kdes[[1]] <- function(x,r=0.00001){
    ct_threshold=40
    y <- numeric(length(x))
    y[x < ct_threshold] <- r
    y[x >= ct_threshold] <- 1-r
    y
  }
  
  den_estimates<-data.frame(matrix(NA,ncol=3,nrow=0))
  colnames(den_estimates)<-c("k","x","y")
  for(i in 1:x){
    z<-replicate(B,sum(sample(viral_loads,i))) #replicate the function B times 
    
    z=z/x
    ct_values<-slope*log10(z)+intercept
   ct_values[ct_values>ct_threshold]=ct_threshold
    #kernsmpooth package:

   kernel_i<-adaptiveKernel(ct_values,to = ct_threshold,from=0)
    all_kdes[[i+1]] <- approxfun(kernel_i$x, kernel_i$y,rule=2)
   den_estimates<-rbind(den_estimates,data.frame("k"= rep(i),"x"=kernel_i$x,"y"=kernel_i$y))
  }
    
  
  return(all_kdes)
}

