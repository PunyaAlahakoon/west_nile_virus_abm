

#load the default theme for figures:
theme<-source("functions/default_theme.R")
default_theme<-default_theme()

source("functions/prev_likelihood_real.R")

#load all the data:
all_data <- readRDS("synthetic_data/all_data.rds")
#load the data: pooled ct values 
cts_per_pool= readRDS("pre_calculations/pooled_cts/cts_per_pool_demo.rds") #they are only data from 2:4 years for the 4 sims, 


ct_threshold=40


#all_dat=subset(all_data,all_data$sim_number==1)
ct_threshold=40

#KDE estimation:
positive_ctss<-all_data %>% filter(ct_value < ct_threshold) %>% pull(ct_value)

#transform to viral loads:
ct_to_vl<-function(ct,intercept,slope){ #ct=m*log10(vl)+c0
  vl<-10^((ct-intercept)/slope)
  vl   
}

vl_to_ct<-function(slope,intercept,vl){ #vl in 10^
  ct1<-slope*log10(vl) +intercept
  #ct<-min(ct1,40)
  ct1
}

intercept=36.9
slope= -2.7 

viral_loads_vec=ct_to_vl(positive_ctss,intercept,slope)

#numbre of pools:

pool_sizes_test=c(1,10, 20, 50,100,200)
f_kde_all=lapply(1:length(pool_sizes_test),function(x) readRDS(file = paste0("f_kde_",pool_sizes_test[x],".rds")))
#read all the fkdes 


n_pools_test=rep(c(1,10,20,50,100,500),each=50)
#n_pools_test=rep(c(50),each=30)
pool_sizes_test=c(1,10, 20, 50,100,200)

#prev=0.5
prev=c(0.001,0.01, 0.1)
n_prev<-length(prev)

couts<-expand.grid(n_pools_test,pool_sizes_test,prev)
colnames(couts)<-c("n_pools_test","pool_sizes_test","prev")

p_tests<-seq(0,1,0.005)
#caclulate pooled ct values and the binarised ddata given the viral loads and the actual prevalances"

ct_pooled_samples=NULL
estimated_prev<-NULL

for(i in 1:nrow(couts)){
  #calculate the kde for each sample size:
  n_positives <-replicate(couts$n_pools_test[i], rbinom(1, couts$pool_sizes_test[i],couts$prev[i]))
  pos_vls<-sapply(1:length(n_positives),function(x) sample(viral_loads_vec,n_positives[x],replace=F))
  neg_vls<-sapply(1:length(n_positives), function(x) couts$pool_sizes_test[i]-n_positives[x])
  pool_vl<-sapply(1:length(n_positives),function(x) mean(c(pos_vls[[x]],rep(0,neg_vls[x]))))
  pool_ct<-vl_to_ct(slope,intercept,pool_vl)
  pool_ct[pool_ct>40]=40 
  ct_pooled_samples=rbind(ct_pooled_samples,data.frame(rep(couts[i,]),"pooled_ct"=pool_ct,"n_positives"=n_positives))
  #estimate the prevalence:
  f_kde=NULL
  idx=which(pool_sizes_test %in% couts$pool_sizes_test[i])
  sub_pools=rep(idx,couts$n_pools_test[i])
  f_kde=f_kde_all[sub_pools]
  
  #  f_kde=NULL
  # f_kde=all_f_kde[rep(couts$pool_sizes_test[i],couts$n_pools_test[i])] #get the kde's sample sizes 
  
  fit <- bbmle::mle2(prev_likelihood_real,start=list(ini_p=runif(1)),
                     data=list(f_kde=f_kde,n_pools=couts$n_pools_test[i],
                               n_sample_size=rep(couts$pool_sizes_test[i],couts$n_pools_test[i]),
                               observed_ct_vec=pool_ct,neg_lik=TRUE),
                     method="Brent",lower=0,upper=1)
  suppressMessages(
    confint1 <- tryCatch({
      ci <- confint(fit)
      if(is.na(ci[1])) ci[1] <- 0
      if(is.na(ci[2])) ci[1] <- 1
      ci
    }, error = function(e) {
      estms<-sapply(1:length(p_tests),function(x) prev_likelihood_real(f_kde=f_kde,ini_p=p_tests[x],n_pools=couts$n_pools_test[i],
                                                                       n_sample_size=rep(couts$pool_sizes_test[i],couts$n_pools_test[i]),
                                                                       observed_ct_vec=pool_ct,neg_lik=FALSE))
      #find the minuim of the negative log likelihood:
      inx=which.max(estms)
      lr<-2*(estms[inx]-estms)
      #calculate the chi -squred related value for the 95% ci:
      chix=max(estms)-(3.84/2)
      #fit <- bbmle::mle2(p_lk,start=list(p_i=0.01),
      # method="Brent",lower=0,upper=1)
      q1<-p_tests[which.min(abs(estms[1:(inx-1)]-chix))]
      q2<-p_tests[(inx+1):length(estms)][which.min(abs(estms[(inx+1):length(estms)]-chix))]
      return(c(q1,q2))
    }, warning = function(w) {
      estms<-sapply(1:length(p_tests),function(x) prev_likelihood_real(f_kde=f_kde,ini_p=p_tests[x],n_pools=couts$n_pools_test[i],
                                                                       n_sample_size=rep(couts$pool_sizes_test[i],couts$n_pools_test[i]),
                                                                       observed_ct_vec=pool_ct,neg_lik=FALSE))
      #find the minuim of the negative log likelihood:
      inx=which.max(estms)
      lr<-2*(estms[inx]-estms)
      #calculate the chi -squred related value for the 95% ci:
      chix=max(estms)-(3.84/2)
      #fit <- bbmle::mle2(p_lk,start=list(p_i=0.01),
      # method="Brent",lower=0,upper=1)
      q1<-p_tests[which.min(abs(estms[1:(inx-1)]-chix))]
      q2<-p_tests[(inx+1):length(estms)][which.min(abs(estms[(inx+1):length(estms)]-chix))]
      return(c(q1,q2))
    })
  )
  
  
  
  #do the CDC estimation as well:
  inf_states<-(pool_ct<40)*1
  pol_sz=rep(couts$pool_sizes_test[i],couts$n_pools_test[i])
  pv<-pooledBin(inf_states,pol_sz)
  
  estimated_prev<-rbind(estimated_prev,data.frame(couts[i,],"prev_pooled"=coef(fit),"pooled_q1"=confint1[1],"pooled_q2"=confint1[2],
                                                  "prev_cdc"=pv$P,"cdc_q1"=pv$Lower,"cdc_q2"=pv$Upper))
}

#saveRDS(estimated_prev, file = "estimated_prev_all_pools_1_10_20.rds")
#saveRDS(estimated_prev, file = "estimated_prev_all_pools_all_4.rds")

estimated_prev_1=readRDS("pre_calculations/estimated_prevalence/estimated_prev_1.rds")
estimated_prev_2=readRDS("pre_calculations/estimated_prevalence/estimated_prev_2.rds")
estimated_prev_3=readRDS("pre_calculations/estimated_prevalence/estimated_prev_3.rds")

estimated_prev=rbind(estimated_prev_1,estimated_prev_2,estimated_prev_3)
estimated_prev=subset(estimated_prev,estimated_prev$n_pools_test==100)
estimated_prev=subset(estimated_prev,estimated_prev$pool_sizes_test %in% c(1,10,20,30,40,50))


m_estimated_prev<-melt(estimated_prev,id=c("n_pools_test","pool_sizes_test","prev"))

dat_fails_cdc=data.frame("pool_sizes_test"=estimated_prev$pool_sizes_test,
                         "prev"=estimated_prev$prev,"prev_cdc"=estimated_prev$prev_cdc)
dat_fails_cdc=dat_fails_cdc[complete.cases(dat_fails_cdc), ]

dat_fails_cx<-dat_fails_cdc %>%
  group_by(pool_sizes_test,prev) %>%
  summarise(mx=max(prev_cdc))
dat_fails_cx

dat_fails_cdx<-dat_fails_cx %>%
  group_by(pool_sizes_test) %>%
  summarise(mxx=max(prev))
dat_fails_cdx

lvs=c("Pool size = 1","Pool size = 10","Pool size = 20","Pool size = 30",
      "Pool size = 40","Pool size = 50")
m_estimated_prev$pool_sizes_test=as.factor(m_estimated_prev$pool_sizes_test)
levels(m_estimated_prev$pool_sizes_test)=lvs

dat_fails_cdx$pool_sizes_test=lvs
dat_fails_cdx$yy=c(NA,0.16,0.15,0.1,0.1,.65)
dat_fails_cdx$xx=c(NA,0.5,0.45,0.35,0.38,0.23)
dat_fails_cdx[1,2]=NA

pl<-ggplot()+
  geom_point(data=m_estimated_prev,aes(x=(prev),y=value,colour =variable),alpha=0.2)+
  scale_color_manual(values = c("#ad7f00","#7c2d1d"),
                     labels=c("Pooled Ct values","Binarised data"))+
  geom_segment(data=dat_fails_cdx,aes(x = mxx,xend = mxx,y=0,yend = 0.7),colour = "brown",linetype = "dotted")+
  geom_label(data=dat_fails_cdx,aes(label=mxx,x = xx,y=yy),colour = "brown",label.size=0.1)+
  facet_wrap(as.factor(n_pools_test)~as.factor(pool_sizes_test))+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.5)+
  default_theme+
  xlim(0,1)+
  xlab("Actual prevalence")+
  ylab("Estimated prevalence")+
  theme(legend.position = "bottom")+
  theme( legend.title = element_blank())

pl  

ggsave("figures/figure_6.png",pl,width = 12,height = 8)

#supplementa figures
#####
#read the files

estimated_prev_1=readRDS("pre_calculations/estimated_prevalence/estimated_prev_all_pools_all_1.rds")
estimated_prev_2=readRDS("pre_calculations/estimated_prevalence/estimated_prev_all_pools_all_2.rds")
estimated_prev_3=readRDS("pre_calculations/estimated_prevalence/estimated_prev_all_pools_all_3.rds")
estimated_prev_4=readRDS("pre_calculations/estimated_prevalence/estimated_prev_all_pools_all_4.rds")

estimated_prev=rbind(estimated_prev_1,estimated_prev_2,estimated_prev_3,estimated_prev_4)

prev_sub<-data.frame("n_pools_test"=estimated_prev$n_pools_test,"pool_sizes_test"=estimated_prev$pool_sizes_test,
                     "prev"=estimated_prev$prev,"prev_pooled"=estimated_prev$prev_pooled,"prev_cdc"=estimated_prev$prev_cdc )

m_prev_sub<-melt(prev_sub,id=c("n_pools_test","pool_sizes_test","prev"))




sot_estimated_prev=subset(prev_sub,prev_sub$n_pools_test!=1 )
sot_estimated_prev=subset(sot_estimated_prev,sot_estimated_prev$prev %in% c(0.001,0.01 ,0.1))

pre_dat=sot_estimated_prev %>%
  group_by(n_pools_test,pool_sizes_test,prev) %>%
  reframe(diff_pooled=prev-prev_pooled,diff_cdc=prev-prev_cdc)

dat<-pre_dat %>%
  group_by(n_pools_test,pool_sizes_test,prev) %>% 
  reframe(mn_pooled=mean(diff_pooled),pooled_q1=quantile(diff_pooled,probs=0.025,na.rm = T),
          pooled_q2=quantile(diff_pooled,probs=0.975,na.rm = T),
          mn_cdc=mean(diff_cdc,na.rm=T),cdc_q1=quantile(diff_cdc,probs=0.025,na.rm = T),
          cdc_q2=quantile(diff_cdc,probs=0.975,na.rm = T) )


lvs=paste0("Pool size = ",c(1,10,20,50,100,200))
dat$pool_sizes_test=as.factor(dat$pool_sizes_test)
levels(dat$pool_sizes_test)=lvs

lv2=paste0("Prevalence = ",c(0.001,0.01,0.1))
dat$prev=as.factor(dat$prev)
levels(dat$prev)=lv2

p1<-ggplot()+
  geom_line(data=subset(dat,prev=="Prevalence = 0.001"),aes(x=as.factor(n_pools_test),y=mn_pooled,group=as.factor(pool_sizes_test)),color="black")+
  geom_ribbon(data=subset(dat,prev=="Prevalence = 0.001"),aes(x=as.factor(n_pools_test),ymin=pooled_q1, 
                                                              ymax=pooled_q2,group=as.factor(pool_sizes_test)),fill="#ee8494", alpha = 0.4)+
  facet_wrap(~as.factor(pool_sizes_test),ncol=1)+default_theme+
  geom_hline(yintercept = 0,linetype="dotted",color="red")+
  ylab("Actual prevalence - estimated prevalence")+
  xlab("Number of pools")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p1


p2<-ggplot()+
  geom_line(data=subset(dat,prev=="Prevalence = 0.01"),aes(x=as.factor(n_pools_test),y=mn_pooled,group=as.factor(pool_sizes_test)),color="black")+
  geom_ribbon(data=subset(dat,prev=="Prevalence = 0.01"),aes(x=as.factor(n_pools_test),ymin=pooled_q1, 
                                                             ymax=pooled_q2,group=as.factor(pool_sizes_test)),fill="#008e7c", alpha = 0.4)+
  facet_wrap(~as.factor(pool_sizes_test),ncol=1)+default_theme+
  geom_hline(yintercept = 0,linetype="dotted",color="red")+
  ylab("Actual prevalence - estimated prevalence")+
  xlab("Number of pools")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p2

p3<-ggplot()+
  geom_line(data=subset(dat,prev=="Prevalence = 0.1"),aes(x=as.factor(n_pools_test),y=mn_pooled,group=as.factor(pool_sizes_test)),color="black")+
  geom_ribbon(data=subset(dat,prev=="Prevalence = 0.1"),aes(x=as.factor(n_pools_test),ymin=pooled_q1, 
                                                            ymax=pooled_q2,group=as.factor(pool_sizes_test)),fill="#816997", alpha = 0.4)+
  facet_wrap(~as.factor(pool_sizes_test),ncol=1)+default_theme+
  geom_hline(yintercept = 0,linetype="dotted",color="red")+
  ylab("Actual prevalence - estimated prevalence")+
  xlab("Number of pools")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p3


p1|p2|p3

ggsave("figures/supp_fig_1.png",last_plot(),height = 10,width = 10)

#subset pooled ct values 
pooled_esti<-data.frame("n_pools_test"=estimated_prev$n_pools_test,"pool_sizes_test"=estimated_prev$pool_sizes_test,
                        "prev"=estimated_prev$prev,
                        "pooled_q1"=estimated_prev$pooled_q1,"pooled_q2"=estimated_prev$pooled_q2)

#check if the cis include the true prevlance 
estimated_prev$within_ct=between(estimated_prev$prev,estimated_prev$pooled_q1,estimated_prev$pooled_q2)*1
estimated_prev$within_cdc=between(estimated_prev$prev,estimated_prev$cdc_q1,estimated_prev$cdc_q2)*1

#get the counts 
sumrry_cts<-estimated_prev %>% 
  group_by(n_pools_test,pool_sizes_test,prev) %>%
  summarise(perce_pooled=(sum(within_ct)/n())*100,perce_cdc=(sum(within_cdc)/n())*100,.groups = 'drop')

m_summary_cts=melt(sumrry_cts,id=c("n_pools_test","pool_sizes_test","prev"))

m_summary_cts=subset(m_summary_cts,m_summary_cts$n_pools_test!=1)
m_summary_cts=subset(m_summary_cts,m_summary_cts$pool_sizes_test!=1)

p2<-ggplot(data=m_summary_cts,aes(x=as.factor(n_pools_test),
                                  y=value,color=as.factor(pool_sizes_test),group = as.factor(pool_sizes_test)))+
  geom_point(size=2,shape=16)+
  facet_wrap(variable~prev)+
  scale_color_viridis(discrete=T,option="D")+
  #scale_colour_manual(
  # values = c("#050ba8", "#5235d5", "#855dff", "#b687ff", "#e6b2ff","#ffe8ff")
  # )+
  default_theme+
  xlab("Number of pools")+
  ylab("Percentage (%)")+
  theme(legend.position = "bottom")+
  labs(color = "Pool size")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

p2

ggsave("figures/pre_coverage.png",last_plot())

p2<-ggplot()+
  geom_line(data=dat,aes(x=as.factor(n_pools_test),y=mn_cdc,group=as.factor(pool_sizes_test),color=as.factor(prev)))+
  geom_ribbon(data=dat,aes(x=as.factor(n_pools_test),ymin=cdc_q1, 
                           ymax=cdc_q2,group=as.factor(pool_sizes_test),fill=as.factor(pool_sizes_test)), alpha = 0.5)+
  facet_wrap(as.factor(pool_sizes_test)~as.factor(prev),ncol=3,scales="free_y")+
  default_theme+
  geom_hline(yintercept = 0,linetype="dotted")+
  theme(legend.position = "none")+
  ylab("Actual prevalence-estimated prevalence")+
  xlab("Number of pools")+
  theme(legend.position = "none")+
  ggtitle(" Estimation with binarised data")

p2



ggsave("figures/pre_test.png",last_plot(),height = 8,width = 14)
