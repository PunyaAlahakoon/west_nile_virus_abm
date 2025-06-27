
#kde approximation for the ct values 
# B number of sims 


#kde function approximation 

##final version: as of 21/10/2024: use the viral loads of the infected mozzies only. use the same transformation used in analysis2
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
  
  #when k=0, the distribution is the distribution to k=1 multiplied by r to account for the false positives 
 # den_estimates<-rbind(den_estimates,data.frame("k"= rep(0),"ct"=0:40,"fk"=rep(0)))
  for(i in 1:x){
    #  sample x viral loads and get the sum 
   # z<-replicate(B,sum(sample(viral_loads,i,replace = F))) #replicate the function B times
    z<-replicate(B,sum(sample(viral_loads,i))) #replicate the function B times 
    
    z=z/x
  # ct_values<-vl_to_ct(-2.7,36.9,z)
    ct_values<-slope*log10(z)+intercept
   ct_values[ct_values>ct_threshold]=ct_threshold
    #kernsmpooth package:

   kernel_i<-adaptiveKernel(ct_values,to = ct_threshold,from=0)
    all_kdes[[i+1]] <- approxfun(kernel_i$x, kernel_i$y,rule=2)
   den_estimates<-rbind(den_estimates,data.frame("k"= rep(i),"x"=kernel_i$x,"y"=kernel_i$y))
  }
    
  #den_estimates$fk[den_estimates$fk<0]=0.01
 #den_estimates[1,3]= den_estimates[1,3]*r
  
  return(all_kdes)
}

