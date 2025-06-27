
#run the em algorithm to estimate the prevalence at each time point
#at the moment, assuming that the sample size is hthe same across all the pools 
#f_kde as a lists of kdes that needs to use pool sizes:\

#this version is for real data, where the pool sizes can be differenet 
prev_likelihood_real<-function(f_kde,ini_p,n_pools,n_sample_size,observed_ct_vec,neg_lik){
 if((0<=ini_p) | (ini_p<=1)){
    
  #p_k|ct,pi-1
  #observed_ct_vec=observed_ct_vec-log2(n_sample_size)
  
  #p_k<-function(fkde,x,observed_ct){
    #f_k<-den_estimates$fk[den_estimates$ct==observed_ct]
   # f_k=rep(NA,length(0:x))
   # for (i in 0:x) {
      #maybe do the following fo reach fkde, give the k value 
    #  fkd_i=subset(fkde,k==i)
     # test<-abs(fkd_i$ct- observed_ct) == min(abs(fkd_i$ct - observed_ct))  #get the closest ct value to the observed value and get the corresponding fk 
      #for each value of x, find the closest observed ct value 
     # f_k[i+1]=fkd_i$fk[test]
    #}
    
    #return(f_k) #return the density given the ct value for all the k values  
  #}
  
  p_lk<-function(p_i,fkdes,b,x_vec,ct_vec){
    pd<-rep(NA,b)
    for (i in 1:b) {
      # pd[i]<-sum((0:x)*p_k(p_i,fkde,x,ct_vec[i]))
      x=x_vec[i]
      fkde=fkdes[[i]]
      pd[i]<-sum(dbinom(0:x,x,p_i)*unlist(lapply(fkde[1:(x+1)], function(f) f(ct_vec[i]))))
                 #  p_k(fkde,x,ct_vec[i]))
    }
    
    lk<-sum(log(pd))
    if(neg_lik==TRUE){
      llk=-(lk) #return the -log likelihood function 
    }else{
      llk=(lk)
    }
    return(llk)
  }
  
  lg_lk<-p_lk(ini_p,f_kde,n_pools,n_sample_size,observed_ct_vec)
  return(lg_lk)}
  
  else{
    lg_lk=0 
  }
}
