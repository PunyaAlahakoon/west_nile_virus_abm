


sample_pools_and_calc_cts<-function(all_data,gridd){
  n_bites=1
#filter the data based on the number of bites/ filter blood seekig mosquitoes. plus just to make sure that you don't have any dead mozzies inclueded as well 
blood_seeking_mozzies<-subset(all_data,all_data$number_of_bites >=n_bites & all_data$inf_state!=2)
#get the population size of the mozzies and then put then into pools 
#prop_sample=.1 #sample size of the mozzies 

#do the sampling based on the observed data:

#read the cscs with pool sizes avg:
avg_sizes<-read.csv("data/sample_sizes_by_surevil_calc.csv")

#get the number of pools per week:
pool_nums<-avg_sizes %>%
  group_by(surv_year,disease_week,state) %>%
  summarise(num_pools=n(),num_pos=sum(WNV))
#pool_nums

#randomly choose a state,given the state get the number of pools, given the state, randomly sample the pool sizes:

pool_sample_dat<-data.frame(matrix(NA,nrow = 0,ncol=5))
colnames(pool_sample_dat)<-c("yr","week","pool_number","pool_size","num_positives")
obs_yrs<-c("2022","2023","2024") 
obs_weeks<-22:39 

obs_grid<-expand.grid(obs_yrs,obs_weeks)

for (i in 1:nrow(obs_grid)) {
  sxij=subset(pool_nums,pool_nums$surv_year==obs_grid[i,1] & pool_nums$disease_week==obs_grid[i,2])
  #x=sample(sxij$state,1) 
  x="Nebraska"
  #get the number of pools for the year and state 
  ni_pools<-sxij$num_pools[sxij$state==x]
  # ni_positives<-sxij$num_pos[sxij$state==x]
  rni_pools<-rpois(1,ni_pools)
  #rni_positives<-sum(rbinom(rni_pools,rni_pools,ni_positives/ni_pools))
  #subset the avg sizes:
  avgij_sizes<-subset(avg_sizes,avg_sizes$surv_year==obs_grid[i,1] & avg_sizes$disease_week==obs_grid[i,2] & avg_sizes$state==x)
  #sample pool sizes for ni_samples:
  p_sizes<-sample(avgij_sizes$pool_size,rni_pools,replace = T)
  pool_sample_dat<-rbind(pool_sample_dat,data.frame("yr"=obs_grid[i,1],"week"=obs_grid[i,2],
                                                    "pool_number"=rni_pools,"pool_size"=p_sizes))
}


#add 0s to the weeks that you don't have
zers_dat<-data.frame("yr"=rep(obs_yrs),"week"=rep(c(1:21,40:52),3),
                     "pool_number"=rep(0),"pool_size"=rep(0))

pool_sample_dat<-rbind(pool_sample_dat,zers_dat)

library(plyr)
pool_sample_dat$yr=revalue(pool_sample_dat$yr,c("2022"="2","2023"="3","2024"=4))
detach("package:plyr", unload = TRUE)



#write a separate function to calculate the pooled ct values:
###add this in a separate r script 
calc_pooled_ct<-function(sim_number,time,week,year,
                         n_pools,r_mozzies,sample_size,mozzies_at_i){
  ct_to_vl<-function(ct,intercept,slope){ #ct=m*log10(vl)+c0
    vl<-10^((ct-intercept)/slope)
    vl   
  }
  
  vl_to_ct<-function(slope,intercept,vl){ #vl in 10^
    ct1<-slope*log10(vl) +intercept
    ct<-min(ct1,40)
    ct
  }
  
  intercept=36.9
  slope= -2.7 
  
  cum_sample<-c(1,cumsum(sample_size))
  
  cts_per_pool=data.frame(matrix(ncol = 6, nrow = 0))
  colnames(cts_per_pool)=c("sim_number","time","week", "year","sample_size", "pooled_ct_per_pool")
  for(j in 1:(n_pools)){
    #generate the sequance of mozzie vector per pool:
    id_seq<-r_mozzies[cum_sample[j]:cum_sample[(j+1)]]
    #get the viral loads of the pool
    pool_vls<-mozzies_at_i$viral_load[mozzies_at_i$id %in% id_seq]
    mn_vl=mean(pool_vls)
    ct_x=vl_to_ct(slope,intercept,mn_vl)
    ct_ss=min(ct_x,40) 
    cts_per_pool<-rbind(cts_per_pool,data.frame("sim_number"=sim_number,"time"=gridd[i,1] ,"week"=week, "year"=year,
                                                "sample_size"=sample_size[j],"pooled_ct_per_pool"=ct_ss))
  }
  cts_per_pool
}


#calculate the pooled ct value for each pool:
cts_per_pool=data.frame(matrix(ncol = 6, nrow = 0))
colnames(cts_per_pool)=c("sim_number","time","week", "year","sample_size", "pooled_ct_per_pool")

#also strore the number of pools separately 
pool_numbers<-data.frame(matrix(ncol = 5,nrow = 0))
colnames(pool_numbers)<-c("sim_number","time","week", "year","number_of_pools")

griddx=subset(gridd,year==yr_to_run & Var2==sim_to_run)
#griddx=gridd

for (i in 1:nrow(griddx)) { #just run one sim for the moment, if not, use 1:nrow(gridd)
  # sample_size=50 #resize the variable all the time 
  
  mozzies_at_i<-subset(blood_seeking_mozzies,blood_seeking_mozzies$current_time==griddx[i,1] & blood_seeking_mozzies$sim_number==griddx[i,2])
  
  #subset pool_sample_dat by week:at the moment, this
  pool_sample_dat_i<-subset(pool_sample_dat,pool_sample_dat$week==griddx$week[i] & pool_sample_dat$yr==griddx$year[i])
  
  #pool sizes per week per year:
  sample_size_at_i<-pool_sample_dat_i$pool_size
  
  n_pools<-pool_sample_dat_i$pool_number[1] #because the pool number is repeated across all the sample sizes 
  if(n_pools>0 ){
    if(length(mozzies_at_i$id)<sum(sample_size_at_i)){#that means the simulated infectious mozzies are not enough to be actually calulated!! 
      cx<-cumsum(sample_size_at_i)
      inx<-which.max(cx[cx<length(mozzies_at_i$id)])
      sample_size_at_i<-sample_size_at_i[1:inx]
      #n_positives=length(positive_samples)
      
    }
    r_mozzies<-sample(mozzies_at_i$id,sum(sample_size_at_i),replace = FALSE)
    n_pools<-length(sample_size_at_i)
    
    cts_pools<-calc_pooled_ct(sim_number=griddx[i,2],time=griddx[i,1],week=griddx[i,3],year=griddx[i,4],
                              n_pools=n_pools,r_mozzies,sample_size_at_i,mozzies_at_i)
    cts_per_pool<-rbind(cts_per_pool,cts_pools)
  }
  pool_numbers<-rbind(pool_numbers,data.frame("sim_number"=griddx[i,2],"time"=griddx[i,1] ,
                                              "week"=griddx[i,3], "year"=griddx[i,4],"number_of_pools"=n_pools))
}

return(cts_per_pool)
}
