#write a function for model for lifespan_function 

lifespan_function<-function(mean_rate,alpha,t){
  t_rate=mean_rate*(1-alpha*sin(pi*t/365))
  return(t_rate)
}
