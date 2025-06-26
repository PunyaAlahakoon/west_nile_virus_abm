

temp_function<-function(mean_rate,alpha,t){
  t_rate=mean_rate*(1-alpha*cos(2*pi*t/365))
  return(t_rate)
}
