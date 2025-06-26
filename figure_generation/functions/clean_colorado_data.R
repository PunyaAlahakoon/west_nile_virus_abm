
clean_colorado_data<-function(){
  col<-read.csv("data/Colorado_Mosquito_Pool_Data_22_23_24.csv")
  
  #col$cq[col$test_code==0]=40
  
  col$test_code[col$cq==55.55]=0
  col$cq[col$cq==55.55]=40
  col$test_code[col$cq<40]=1
  
  
  dat_col<-data.frame("surv_year"=as.factor(col$year),"disease_week"=col$week,"species"=col$spp,
                      "pool_size"=col$total,"WNV"=col$test_code, "ctval"=col$cq)
  dat_col$state<-rep("Colorado")
  dat_col=subset(dat_col,dat_col$surv_year!=2024)
  
  return(dat_col)
}