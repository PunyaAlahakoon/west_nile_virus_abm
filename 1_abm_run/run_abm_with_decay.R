#the following script runs the abm. 
#change parameters as needed in this script. But this code with the parameters will generate synthetic data that contain in the synthetic_data folder. 
#do not change any of the other files in this folder 

#this version of the model assumes productive infection:
#1. 100% of the virus is trasffered to the mozzies after biting a bir
#2. there is 50-50 chance that the mozzie virus dynamics will end up with non-productive and productive infection 

sourceCpp("interface.cpp") #load the interface file, this file has the ABM in C++ 

#number of years to run the sim 

yrs<-4
time=365*yrs#days 
nHosts=1000
nMozzies=10000 
probHostInfection=0.88
#probMozzInfection=0.5
iniInfectedMozzies=100
mozBirthRate=(0.537*10000)*0.054 
#calculations for the birth rate: 0.537 is the per capita birth rate (assume this is the rate of eggs being layed). 
#assume 0.537*10000 is the rate of eggs being layed per day.
#0.054 is the proportion of mozzies surviving from eggs to adults . so the overall birth rate is (0.537*10000)*0.054
egg_to_adult_days=8 #(7-10 days) from cdc 
#https://www.cdc.gov/mosquitoes/about/life-cycle-of-culex-mosquitoes.html#:~:text=Culex%20species%20mosquitoes,develop%20into%20an%20adult%20mosquito. 


host_birth_rate=0.023*nHosts 
host_death_rate=0.0015 #birds 
vertical_transms_prob=0 #when the births include vertical transmission. at the moment, don't include vertical transmission. assume zero 


seed_mozzies=FALSE #if true, then add additional vertical transmission at each time step including the vertical tranmssion prob at birth 
seed_inf_mozzie_rate=0 #seed mozzies to the population, no seeding at the moment 
#have an on-off step for seeding with time: input as a vector 
seed_on<-rep(c(rep(0,150),rep(1,30),rep(0,366-180)),yrs)*seed_inf_mozzie_rate 


source("functions/temp_function.R")
source("functions/lifespan_function.R")
#mozzies: 
mozzie_birth_rate_t=temp_function(mean_rate=mozBirthRate,alpha=1,t=0:time+1)
probMozzInfection_t=temp_function(mean_rate = 0.19,alpha=0.79,t=0:time+1)#at the moment this is only for one species of mozzies
daily_biting_rate=temp_function(mean_rate = .17,alpha=0.17,t=0:time+1)
#plot(daily_biting_rate)
#daily_life_span=temp_function(mean_rate = 100,alpha=0.9,t=0:time+1)
#daily_biting_rate=temp_function(mean_rate = .1,alpha=.5,t=0:time+1)
dt=1 #time step 
daily_pr_bite=1-exp(-daily_biting_rate*dt)
#plot(daily_pr_bite)

#assume that the probability of death is also temperature dependent 
rate_death=temp_function(mean_rate =0.0145, alpha=1,t=0:time+1) #maximum of the death rate at high temperature is approximately 0.029 
#plot(rate_death)
daily_pr_death=1-exp(-rate_death)
#plot(daily_pr_death)
#daily_pr_death=rep(0.029,731) #for constant death rate 


mozzie_viral_dynamics_method<-2
#if 1= mozzies follow their own ct model despite the birds)
#if 2= mozzies use the bird ct value and remain the same throughout the time OR
# mozzies use the ct value of the bird as the min ct and increase from there 

#overwintering period
set_the_overwintering_period<-rep(c(rep(1,63),rep(0,(181-30)),rep(1,(122+30))),yrs)
set_the_mozzie_active_period<-rep(c(rep(0,63),rep(1,(181-15)),rep(0,(122+15))),yrs)

mozzie_owerwintering_prob<-rep(0.25*(1-1*sin(pi*(0:365)/365)),yrs)

active_prob<-0.05*(1-1*cos(pi*(0:(181-15)*2)/(181-15)))


mozzie_becoming_active_prob<-rep(c(rep(0,63),active_prob,rep(0,(122+15))),yrs)




#PART X: use these parameters if you want to test different viral inheritance combinations 
# change_prob_model=c(0,.01,0.1,.5,.99,1) # probability of non-productive infection 
# vira_per_ch<-c(1,.99,.5,.1,0.01,0) # probability of the amount of virus transferred to the mozzie from the bird 

change_prob_model=0.5 # probability of non-productive infection 
vira_per_ch<-1 # probability of the amount of virus transferred to the mozzie from the bird 

#add another factor for decays;
# high_decay_rate<-0.9
# medium_decay_rate<-0.5
# low_decay_rate<-0.1
# decay_rates<-c(0,low_decay_rate,medium_decay_rate,high_decay_rate)
decay_rates<-c(0,0.00001,0.001,0.01,0.05,0.1,0.5,0.9)
grds<-expand.grid(change_prob_model,vira_per_ch,decay_rates)

n_sims<-nrow(grds)#number of sims to run 

for (i in 1:n_sims) { #if running PART X:, chnage n_sims to nrow(grds)
  #if running PART X: 
  #prob_ct_model_change<-grds[i,1]#if the mozzie_viral_dynamics_method==1
  #viral_percent_from_birds<-grds[i,2]
  
  prob_ct_model_change<-grds[i,1]#if the mozzie_viral_dynamics_method==2
  viral_percent_from_birds<-grds[i,2]
  decay_rate<-grds[i,3]
  
  #run the sim several times:
  run_mozzie_model( time,  nHosts,  nMozzies,  
                    daily_pr_bite,
                    mozzie_birth_rate_t,
                    probHostInfection,  
                    probMozzInfection_t, 
                    iniInfectedMozzies,
                    host_birth_rate,
                    host_death_rate,
                    vertical_transms_prob,
                    daily_pr_death,
                    seed_on,seed_mozzies,
                    i,
                    egg_to_adult_days,
                    decay_rate,
                    mozzie_viral_dynamics_method, 
                    prob_ct_model_change,
                    viral_percent_from_birds,
                    set_the_overwintering_period,
                    mozzie_owerwintering_prob,
                    set_the_mozzie_active_period,
                    mozzie_becoming_active_prob,
                    FALSE)
  print(memuse::Sys.meminfo())
  gc()
}


