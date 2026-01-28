

library(car)
library(ggplot2)
library(reshape2)
library(cowplot)
library(ggpubr)
library(ggdist)
library(viridis)
library(viridisLite)
library(colorspace)
library(wesanderson)
library(PooledInfRate)
library(dplyr)
library(HDInterval)
library(foreach)
library(doParallel)
library(mixtools)



default_theme<-function(){
  
  theme=theme_minimal(15)+
    theme(plot.background=element_blank(),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          #axis.title.y=element_blank(),
          #strip.text.x = element_text(size = 14),
          # strip.text.y = element_text(size = 14)
    )+
    theme(panel.grid.major.y = element_line(color = "grey",
                                            linewidth = 0.5,
                                            linetype = "dotted"))+
    theme(legend.title=element_blank(),legend.background = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(theme)
  
}

#load the data;
# read the weekly data
yrs<-4
time=365*yrs#days 

sample_size=50
#max_pools=100

n_bites=0 #number of bites to consider 


d_lod=40  #limit of detection of the ct value 



fq=7#how frequent the ct values to be observed 
sample_frq=seq(fq,time,by=fq) #start from the 7th day 
#remove the first year ,comment this if you don't need this:
sample_frq=sample_frq[-(1:52)]


change_prob_model=0.5 # probability of non-productive infection 
vira_per_ch<-1 # probability of the amount of virus transferred to the mozzie from the bird 

# #add another factor for decays;
# high_decay_rate<-0.9
# medium_decay_rate<-0.5
# low_decay_rate<-0.1
# decay_rates<-c(0,low_decay_rate,medium_decay_rate,high_decay_rate)
decay_rates<-c(0,0.00001,0.001,0.01,0.05,0.1,0.5)
#decay_rates<-0
grds<-expand.grid(change_prob_model,vira_per_ch,decay_rates)

grds<-expand.grid(change_prob_model,vira_per_ch,decay_rates)

n_sims<-nrow(grds)#number of sims to run 

days<-rep(sample_frq,each=n_sims)
grds$n_sims<-1:n_sims

gridd<-data.frame(days,rep(grds[,4]),rep(grds[,1]),rep(grds[,2]),rep(grds[,3])) #repeats each of the decay rates vector 
gridd$week<-rep(rep(1:52,(yrs-1)),each=n_sims) #remove the first year 
gridd$year<-rep(rep(2:yrs,each=52),each=n_sims)
colnames(gridd)<-c("days","n_sims","change_prob_model","vira_per_ch","decay_rates","week","year")

#gridd<-subset(gridd,decay_rates %in% c(0,0.1))


pop_size=c()
all_data<- vector(mode="list",length=length(sample_frq))

for(i in 1:nrow(gridd)){
  # mozzie_objects_" << day << "_"<< sim_number <<"_" << prob_ct_model_change << "_" << viral_percent_from_birds << "_" << decay_rate << ".csv";

  filename<-paste0("mozie_outputs_wth_decay_model/mozzie_object_data/mozzie_objects_", gridd[i,1], "_", 
                   gridd[i,2], "_", gridd[i,3], "_", gridd[i,4], "_", gridd[i,5], ".csv")
  dat=read.csv(filename)
  dat=subset(dat,dat$inf_state!=2)
  dat$week=gridd[i,6]
  dat$year=gridd[i,7]
  dat$prob_model<-gridd[i,3]
  dat$viral_percent<- gridd[i,4]
  dat$decay_rate<-gridd[i,5]
  
  all_data[[i+1]] <- dat
  pop_size=c(pop_size,nrow(dat))
  
}

all_data <- do.call("rbind",all_data)

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

#blood_seeking_mozzies<-subset(all_data,all_data$number_of_bites >=n_bites & all_data$inf_state!=2)

#
pools<-sample_pools_and_calc_cts(all_data,gridd)

all_pools_lees_40<-subset(pools,pools$pooled_ct_per_pool<40)

#get pools that are less than 40 


cts_per_pool= readRDS("data_synthetic/cts_per_pool.rds") 
deca<-0
no_dec_cts_per_pool<-subset(cts_per_pool,cts_per_pool$pooled_ct_per_pool<40)
no_dec_cts_per_pool<-cbind("decay_rate"=rep(deca),no_dec_cts_per_pool[,2:6])


all_pools_lees_40<-subset(all_pools_lees_40,decay_rate!=0)
all_pools_lees_40<-rbind(no_dec_cts_per_pool,all_pools_lees_40)


all_pools_lees_40$decay_rate=as.factor(all_pools_lees_40$decay_rate )
levels(all_pools_lees_40$decay_rate)<-paste0("Decay rate = ",decay_rates)

all_pools_lees_40$year=as.factor(all_pools_lees_40$year)
levels(all_pools_lees_40$year)<-c("Year 1","Year 2","Year 3")

p1 <- ggplot(data = all_pools_lees_40,
             aes(x = pooled_ct_per_pool, fill = as.factor(year))) +
  geom_density(adjust=2,alpha = 0.6) +
  facet_wrap(~ as.factor(decay_rate), ncol = 2,scales="free") +
  xlim(0, 40) +
  scale_fill_manual(
    values = c("#0085aa", "#816997", "#d69b56", "#00847f", "#ee8494")
  ) +
  labs(fill = "Year") +          # <-- Change legend title here
  default_theme()+
  xlab("Pooled Ct values \n (less than the LOD)") +
  ylab("Density") +
  theme(legend.position = "bottom")
p1

ggsave("figures/figure_decay_rates.png",p1, width = 10, height = 10)


p_ct_dynamics+p1

