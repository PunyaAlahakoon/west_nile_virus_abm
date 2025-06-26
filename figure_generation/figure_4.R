
#load the default theme for figures:
theme<-source("functions/default_theme.R")
default_theme<-default_theme()

# read the weekly data
yrs<-4
time=365*yrs#days 
n_sims<-4#number of sims 

sample_size=50
#max_pools=100

n_bites=0 #number of bites to consider 

#out of the sims that's been run, which one to visulaise and use for prevalance estimation? 
yr_to_run=4
sim_to_run=2

#KDE estimation related parameters 
epsilon=0.0005
bandwidth=0.1
mean_r=c(40,42)
std<-c(0.1,1)
d_lod=40  #limit of detection of the ct value 



fq=7#how frequent the ct values to be observed 
#sample_frq=seq(fq,time,by=fq) #start from the 7th day 
#remove the first year ,comment this if you don't need this:
#rem_y<-364*(2:yrs)+7
sample_frq=c()
fq=7
s_frq=seq(fq,365,by=fq)
sample_frq=c(sample_frq,s_frq)

for (i in 2:yrs) {
  s_frq=s_frq+365
  sample_frq=c(sample_frq,s_frq)
}
sample_frq=sample_frq[-c(1:52)]

#sample_frq=c(seq(fq,365,by=fq),  )
#
#sample_frq=sample_frq[!sample_frq %in% rem_y]

#read the mazzie states and hosts and loop through them 
mozzie_states=vector(mode="list",length=length(n_sims))
host_states=vector(mode="list",length=length(n_sims))
for(i in 1:n_sims){
  #mozzi_filename=paste0("mozzie_outputs/states/mozzie_states_",i,"_0.99_1.csv")
  mozzi_filename=paste0("~/Documents/GitHub/west_nile_virus_abm/abm_run/mozzie_outputs/states/mozzie_states_",i,"_0.5_1.csv")
  dat1<-read.csv(mozzi_filename)
  dat1$sim_number=i
  dat1=dat1[-1,]
  dat1$time=1:(365*yrs)
  dat1$year=rep(1:yrs,each=365)
  dat1=subset(dat1,time %in% sample_frq)
  #dat1$sim_number=1
  mozzie_states[[i+1]]<-dat1
  #host_filename=paste0("host_outputs/states/host_states_",i,"_0.99_1.csv")
  host_filename=paste0("~/Documents/GitHub/west_nile_virus_abm/abm_run/host_outputs/states/host_states_",i,"_0.5_1.csv")
  dat2<-read.csv(host_filename)
  dat2$sim_number=i
  dat2=dat2[-1,]
  dat2$time=1:(365*yrs)
  dat2$year=rep(1:yrs,each=365)
  dat2=subset(dat2,time %in% sample_frq)
  
  #dat2$week=rep(rep(1:52,each=7),yrs)
  #dat2$sim_number=1
  host_states[[i+1]]<-dat2
}
mozzie_states<-do.call("rbind",mozzie_states)
host_states<-do.call("rbind",host_states)

mozzie_states=mozzie_states[,-8] 
host_states=host_states[,-8] 


#add weeks and year columns 
host_states$weeks<-rep(1:52,length.out=nrow(host_states))
#host_states$year<-rep(2:yrs,each=52)
mozzie_states$weeks<-rep(1:52)
mozzie_states$year<-rep(2:yrs,each=52)

#all_data=data.frame()
#load the mozzie object data with values 
pop_size=c()
all_data <- vector(mode="list",length=length(sample_frq))
gridd<-expand.grid(sample_frq,1:n_sims) #run this when running assumption 1
#gridd<-data.frame(sample_frq,rep(1)) #this is just the sim number under assumption 2

gridd$week<-rep(rep(1:52,(yrs-1)),n_sims) #remove the first year 
gridd$year<-rep(rep(2:yrs,each=52),n_sims)

for(i in 1:nrow(gridd)){
  filename<-paste0("~/Documents/GitHub/west_nile_virus_abm/abm_run/mozzie_outputs/mozzie_object_data/mozzie_objects_", gridd[i,1], "_", gridd[i,2], "_0.5_1.csv")
  dat=read.csv(filename)
  dat=subset(dat,dat$inf_state!=2)
  dat$week=gridd[i,3]
  dat$year=gridd[i,4]
  all_data[[i+1]] <- dat
  pop_size=c(pop_size,nrow(dat))
}

all_data <- do.call("rbind",all_data)
saveRDS(all_data, "all_data.rds")

#mozzie_states[,1:5]=mozzie_states[,1:5]
mozzie_states[,1:5]=mozzie_states[,1:5]/mozzie_states[,6]
#mozzie_states=mozzie_states[,-6] #remove the pop because it doesn't really matter anymore because the prevalance is per current population 
mozzie_states=mozzie_states[,-5] #remove the new adults  

host_states[,1:5]=host_states[,1:5]/host_states[,6]
#host_states=host_states[,-6] #remove the pop because it doesn't really matter anymore because the prevalance is per current population 


names(mozzie_states)[names(mozzie_states)=="Total_pop"]<-"Population size"

m_host_states<-melt(host_states,id=c("time","sim_number","weeks","year"))
m_mozzie_states<-melt(mozzie_states,id=c("time","sim_number","weeks","year"))

#remove the first year of the data:
yr_to_keep=365:(365*yrs)
####to edit 
m_host_states=subset(m_host_states,m_host_states$time  %in% yr_to_keep)
m_mozzie_states=subset(m_mozzie_states,m_mozzie_states$time %in% yr_to_keep)

m_host_states$sim_id<-1:nrow(m_host_states)


p_hosts_2<-ggplot(data=m_host_states,aes(x=weeks,y=value,group=interaction(sim_number,year)))+
  facet_wrap(~variable,ncol = 3,scales = "free")+
  #  geom_point(size=0.5)+
  geom_line(color="#2f2113",alpha=0.5)+
  default_theme+
  #  scale_x_continuous(breaks=seq(7, 1820, 28*6))+ 
  xlab("Time (days)")+
  ylab("Prevalance \n (per current population)")+
  #scale_color_manual(values = "#816997")+
  #scale_color_manual(values  = wes_palette("BottleRocket1", 3900, type = "continuous"))+
  theme(legend.position = "none")
p_hosts_2


ggsave("figures/host_states_3_yrs_5_sims.png",p_hosts_2,width = 8,height = 5)

levels(m_mozzie_states$variable)=c("Susceptible","Infectious","Deaths","Births","Population_size","Daily_incidence")

m_mozzie_states$variable<-factor(m_mozzie_states$variable,levels = c("Susceptible","Infectious","Daily_incidence","Population_size",
                                                                     "Deaths","Births") )

m_mozzie_states$sim_id=1:nrow(m_mozzie_states)

# 

p_mozzies_2<-ggplot(data = m_mozzie_states,aes(x=weeks,y=value,group=interaction(sim_number,year),color=variable))+
  facet_wrap(~variable,ncol = 3,scales = "free_y",strip.position = "left", 
             labeller=as_labeller(c(Susceptible="Susceptible",Infectious="Prevalence",
                                    Deaths="Daily deaths",Births="Daily births",Population_size="Population size",Daily_incidence="Daily incidence")))+
  geom_line(alpha=0.5,)+
  # geom_point(size=0.5)+
  #  geom_hline(yintercept=0,linetype="dashed")+
  # geom_hline(yintercept=0.05,linetype="dashed")+
  # geom_vline(xintercept = c(24,39))+
  scale_color_manual(values=c("#c2620d","#c3454b","#a44372","#714c81","#414f75","#2f4858"),labels=c("Susceptible","Infectious","Daily_incidence","Population_size",
                                                                                                    "Deaths","Births") )+
  #  scale_color_manual(values = "#005b48")+
  xlab("Time (weeks)")+
  ylab(NULL)+
  # ylab("Prevalance \n (per current population)") +
  default_theme+
  theme(strip.placement = "outside")+
  theme(legend.position = "none")

p_mozzies_2

ggsave("figures/mozzie_states_3_yrs_5sims.png",p_mozzies_2,width = 8,height = 5)

summry<-mozzie_states %>% 
  group_by(sim_number,year) %>%
  summarise(max_prev=max(Infectious),max_week=weeks[which.max(Infectious)],max_inc=(max(new_inf_mozzies)/Susceptible[which.max(new_inf_mozzies)])*100)
summry


inf_moz<-all_data$inf_state==1

all_data$time_since_infection=NA
all_data$time_since_infection[inf_moz]=all_data$current_time[inf_moz]-all_data$infected_time[inf_moz]


#plot the following. based one one sim only;

fig_4_manscpt_dat<-subset(all_data,all_data$year==yr_to_run & all_data$sim_number==sim_to_run)


p_mozzie_prev<-ggplot()+
  geom_line(data=subset(m_mozzie_states,variable=="Infectious" & time!=0 & m_mozzie_states$year!=yr_to_run & m_mozzie_states$sim_number!=sim_to_run),
            aes(x=(weeks),y=value,group=interaction(sim_number,year)),color="gray")+
  geom_line(data=subset(m_mozzie_states,variable=="Infectious" & time!=0 & m_mozzie_states$year==yr_to_run & m_mozzie_states$sim_number==sim_to_run)
            ,aes(x=(weeks),y=value,group=1),color="#2f2113",size=1)+
  xlab("Time (weeks)")+
  #geom_vline(xintercept=30,linetype="dotdash")+
  #geom_vline(xintercept=11,linetype="dotdash",color="blue")+
  #geom_vline(xintercept=13,linetype="dotdash")+
  # geom_hline(yintercept=0.05,linetype="dashed")+
  ylab("Prevalence \n (per current population)") +
  ggtitle("(A) Mosquito infection prevalence")+
  default_theme 
p_mozzie_prev

p_inci<-ggplot()+
  geom_line(data=subset(m_mozzie_states,variable=="Daily_incidence" & time!=0 & m_mozzie_states$year!=yr_to_run & m_mozzie_states$sim_number!=sim_to_run),
            aes(x=(weeks),y=value,group=interaction(sim_number,year)),color="gray")+
  geom_line(data=subset(m_mozzie_states,variable=="Daily_incidence" & time!=0 & m_mozzie_states$year==yr_to_run & m_mozzie_states$sim_number==sim_to_run)
            ,aes(x=(weeks),y=value,group=1),color="#2f2113",size=1)+
  xlab("Time (weeks)")+
  #geom_vline(xintercept=30,linetype="dotdash")+
  #geom_vline(xintercept=11,linetype="dotdash",color="blue")+
  #geom_vline(xintercept=13,linetype="dotdash")+
  # geom_hline(yintercept=0.05,linetype="dashed")+
  ylab("Incidence") +
  ggtitle("(B) Mosquito daily incidence")+
  default_theme 
p_inci

tsi_prev<-p_mozzie_prev/p_inci
tsi_prev

pct_l<-ggplot(data=fig_4_manscpt_dat,aes(x=week,y=ct_value))+
  geom_point(color="#2f2113",alpha=0.2,size=0.8)+
  scale_y_reverse()+
  default_theme+
  ylim(40,0)+
  xlab("Time (weeks)")+
  ylab("Ct values of the \n mosquitoes in the population")+
  ggtitle("(C) Individual Ct values")
pct_l

#ggsave("figures/all_cts_model_1.png",pct_l)

#subset the mozzies based on their characteristics

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



all_data$viral_load<-ct_to_vl(all_data$ct_value,intercept,slope )

hist(log10(all_data$viral_load[all_data$ct_value<40]))

hist((all_data$ct_value[all_data$ct_value<40]))




xx=subset(all_data,all_data$ct_value<40)

dx<-xx %>%
  group_by(ct_value) %>%
  filter(n() >1) %>%
  ungroup()
dx

yy=unique(xx$ct_value)
rpeated_ids=xx[!yy,]

#all_data$viral_load[all_data$inf_state==1 & all_data$ct_value==40]=10^2
all_data$viral_load[all_data$inf_state==0 & all_data$ct_value==40]=0



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
pool_nums


#pool_avg<-pool_nums %>%
# group_by(disease_week) %>%
# summarise(mean_num_pools=round(mean(num_pools)))
#pool_avg

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


saveRDS(cts_per_pool, "cts_per_pool.rds")


ct_per_pools<-ggplot(data=subset(cts_per_pool,cts_per_pool$year==yr_to_run & cts_per_pool$sim_number==sim_to_run),aes(x=(week),y=pooled_ct_per_pool))+
  #scale_y_reverse()+
  geom_point(colour="#bd871f",alpha=0.6)+
  ylim(40,0)+
  xlim(0,52)+
  scale_fill_manual(values  = wes_palette("AsteroidCity1", 3, type = "continuous"))+
  default_theme+
  theme(legend.position = "bottom")+
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  xlab("Time (weeks)")+
  ylab("Calculated Ct \n values per pool")+
  ggtitle("(D) Pooled Ct values")

ct_per_pools


tsi_prev_ct<-(p_mozzie_prev + p_inci)/(pct_l+ct_per_pools)
tsi_prev_ct

ggsave("figures/figure_4.png",tsi_prev_ct, width = 12, height = 8)




