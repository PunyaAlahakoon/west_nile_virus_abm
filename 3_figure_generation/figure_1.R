#Figure 1:

#load the default theme for figures:
theme<-source("functions/default_theme.R")
default_theme<-default_theme()


#load the data:
nebarska<-source("functions/clean_nebraska_data.R")
dat_neb<-clean_nebraska_data()

#load colorado data:
colorado<-source("functions/clean_colorado_data.R")
dat_col<-clean_colorado_data()


#conbine both datasets:
all_dat<-rbind(dat_neb,dat_col)

lvs<-c("Nebraska","Colorado") #change the order of the states when visualising 
all_dat$state=factor(all_dat$state,levels = lvs)

all_dat$pool_size[all_dat$pool_size>50]=50

sum<-subset(all_dat,state=="Colorado") %>%
  group_by(disease_week) %>%
  reframe(counts=n(),mean_ct=mean(ctval,na.rm=T)
          ,min_ct=min(ctval,na.rm=T),max_ct=max(ctval,na.rm=T),
          n_postives=sum(WNV),men_ct_less_40=mean(subset(ctval,ctval<40),na.rm=t),
          avg_pool_sz=ceiling(mean(pool_size)))
sum


sum<-subset(all_dat,pool_size<=25 & ctval!=40) %>%
  group_by(state) %>%
  reframe(counts=n(),mean_ct=mean(ctval,na.rm=T),min_ct=min(ctval,na.rm=T),max_ct=max(ctval,na.rm=T))
sum

sumry=all_dat %>% 
  group_by(state) %>%
  summarise(mean=mean(ctval),min=min(ctval),max=max(ctval))
sumry

summary(all_dat)
p_v1<-ggplot(all_dat,aes(x=disease_week,y=ctval,colour=state))+
  geom_point()+
  scale_y_reverse()+
  facet_wrap(as.factor(state)~(as.factor(surv_year)),ncol=3,scales="free_x")+
  default_theme+
  xlim(20,40)+
  xlab("Disease week")+
  ylab("Pooled Ct values")+
  scale_color_manual(values  =c("#9d88b2","#a58000"))+
  # scale_color_manual(values  =c("#9d88b2","#a58000"))+
  theme(legend.position = "none")
#  theme(plot.margin = margin(.5,.5,.5,.5,"cm"))

p_v1



summary_dat<-all_dat %>% 
  group_by(disease_week,state,surv_year) %>%
  summarise(pcr_positive=sum(WNV),pcr_negative=sum(WNV==0))
summary_dat



m_summary_dat<-melt(summary_dat,id=c("disease_week","state","surv_year"))

#plot the binarised data:
p_bin<-ggplot(data=m_summary_dat,aes(x=(disease_week),y=value,fill=variable))+
  geom_bar(stat="identity")+
  xlab(as.numeric(m_summary_dat$disease_week))+
  scale_x_continuous(limits = c(20,40))+
  #scale_x_discrete(breaks=as.factor(c(0,20,25,30,35,40)))+
  # scale_x_discrete(labels=abbreviate)+
  #geom_point()+
  #  geom_count()+
  facet_wrap(as.factor(state)~(as.factor(surv_year)),ncol=3,scales="free_x")+
  default_theme+
  xlab("Disease week")+
  ylab("Pooled RT-qPCR tests")+
  scale_fill_manual(values  =c("#903920","#759228"),labels= c("PCR positive", "PCR negative"))+
  theme(legend.position = "bottom")

p_bin


#load data from human cases:
human_cases<-read.csv("data/human_cases_99_24.csv")
lvs<-c("Jan",
       "Feb",
       "Mar",
       "Apr",
       "May",
       "Jun",
       "Jul",
       "Aug",
       "Sep",
       "Oct",
       "Nov",
       "Dec")
human_cases$Month=factor(human_cases$Month,levels = lvs)

p_human<-ggplot(data=human_cases,aes(as.factor(Month),y=Reported.Cases))+
  geom_bar(stat="identity",fill="#386da5")+
  geom_text(aes(label=Reported.Cases), vjust=0.0001, color="black",
            position = position_dodge(0.9), size=3.5)+
  default_theme+
  xlab("Month of illness onset")+
  ylab("Reported cases")+
  ggtitle("(A) WNV human cases from 1999 - 2024")+
  theme(plot.margin = margin(.5,.5,.5,.5,"cm"))

p_human

# Mosquito 
p_obs<-(p_v1+ggtitle("(A) Mosquito data: Ct values from 2022 - 2024"))/(p_bin+ggtitle("(B) Mosquito data: RT-qPCR test results from 2022 - 2024"))
p_obs

#p_all_obs<-(p_human|plot_spacer())/(p_obs) +plot_layout(heights = c(1, 4))
#p_all_obs

#ggsave("figures/p_all_obs_human.png",p_all_obs,width = 12,height = 15)

ggsave("figures/figure_1.png",p_obs,width = 10,height = 14)

