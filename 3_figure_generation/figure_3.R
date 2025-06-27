

#load the default theme for figures:
source("functions/default_theme.R")
default_theme<-default_theme()


cts_per_pool= readRDS("pre_calculations/pooled_cts/cts_per_pool_assumption_1.rds")
cts_per_pool$sample_size[cts_per_pool$sample_size>50]=50
cts_per_pool_ass1<-subset(cts_per_pool,cts_per_pool$sim_number==25)
cts_per_pool_ass1$assumption=rep("(A) Simulated data: \n Mosquitoes undergo \n dynamic infection from Ct=40")

cts_per_pool= readRDS("pre_calculations/pooled_cts/cts_per_pool.rds") #they are only data from 2:4 years for the 4 sims, 
cts_per_pool_ass4=subset(cts_per_pool,cts_per_pool$sim_number==3)
#also doubke check if the sample size is greater than 50:
cts_per_pool_ass4$sample_size[cts_per_pool_ass4$sample_size>50]=50

cts_per_pool_ass4$assumption=rep("(B) Simulated data: \n Mosquitoes inherit 100% \n of viral load from birds; \n 50% of mosquitoes undergo a \n static viral kinetics trajectory")

all_cts_ass<-rbind(cts_per_pool_ass1,cts_per_pool_ass4)
all_cts_ass_40=subset(all_cts_ass,all_cts_ass$pooled_ct_per_pool<40)

#load the data:
source("functions/clean_nebraska_data.R")
dat_neb<-clean_nebraska_data()

neb_dat<-data.frame("year"=dat_neb$surv_year,"week"=dat_neb$disease_week,"pooled_ct"=dat_neb$ctval,"assumption"=rep("(C) Nebraska data"))
neb_dat$year=as.factor(neb_dat$year)
levels(neb_dat$year)=c("2","3","4")


#also add colorado data, because nebraska doesn't seem to be enough 
#load colorado data:
source("functions/clean_colorado_data.R")
dat_col<-clean_colorado_data()


col_dat<-data.frame("year"=dat_col$surv_year,"week"=dat_col$disease_week,"pooled_ct"=dat_col$ctval,"assumption"=rep("(D) Colorado data"))
col_dat$year=as.factor(col_dat$year)
levels(col_dat$year)=c("2","3","4")

cts_per_pool_dat<-rbind(neb_dat,col_dat)
cts_per_pool_assums_40_dat<-subset(cts_per_pool_dat,cts_per_pool_dat$pooled_ct<40)
levels(cts_per_pool_assums_40_dat$year)<-c("Year 1","Year 2","Year 3")


all_cts_ass_dat<-data.frame(all_cts_ass[,c(4,3,6,7)])
colnames(all_cts_ass_dat)<-c("year","week","pooled_ct","assumption")
all_cts_ass_2<-rbind(all_cts_ass_dat,cts_per_pool_dat)
all_cts_ass_40=subset(all_cts_ass_2,all_cts_ass_2$pooled_ct<40)
all_cts_ass_40$year=as.factor(all_cts_ass_40$year)
levels(all_cts_ass_40$year)<-c("Year 1","Year 2","Year 3")

p3<-ggplot(data=all_cts_ass_40,aes(x=pooled_ct,fill = as.factor(year)))+
  geom_density(alpha=0.6)+
  facet_wrap(~assumption,ncol=2,scales="free_x")+
  xlim(0,40)+
  ylim(0,.1)+
  scale_fill_manual(values  = c("#0085aa","#816997","#d69b56","#00847f","#ee8494"))+
  default_theme+
  # ylim(0,40)+
  xlab("Pooled Ct values \n (less than the LOD)")+
  ylab("Density")+
  theme(legend.position = "bottom")
p3


ggsave("figures/figure_3.png",last_plot(),height = 8,width = 10)
