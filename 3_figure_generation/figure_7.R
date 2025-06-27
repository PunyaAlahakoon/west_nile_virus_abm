
#load the default theme for figures:
source("functions/default_theme.R")
default_theme<-default_theme()

#read all the fkdes 
all_f_kde=NULL
all_f_kde=sapply(1:50,function(x) readRDS(file = paste0("pre_calculations/kde_calcs/f_kde_",x,".rds")))

#read the clean real datasets 
source("functions/clean_nebraska_data.R")
dat_neb<-clean_nebraska_data()
source("functions/clean_colorado_data.R")
dat_col<-clean_colorado_data()

source("functions/prev_likelihood_real.R")


###additional estimates 
p_tests<-seq(0,.5,0.005)

#extract the disease weeks and years 
yrs<-2022:2024
weeks<-24:39

#repeat the procedure for each year:
prev_esti_by_ys_neb<-NULL

#there's a point that has more than 50 num_count:
dat_neb$pool_size[dat_neb$pool_size>50]=50


for (j in 1:3) {
  yr_to_run=yrs[j]
  cts_of_yr<-subset(dat_neb,surv_year==yr_to_run)
  #run the prevalence estimation in parallel:
  optim_prev_esti<-NULL
  for(i in 1:52) {
    cts_at_i<-subset(cts_of_yr,cts_of_yr$disease_week==i)
    #remove the row with NA ct values 
    cts_at_i=cts_at_i[complete.cases(cts_at_i),]
    b_pools=nrow(cts_at_i)
    if(b_pools>0){
      sample_size_i=cts_at_i$pool_size
      observed_ct_vec=cts_at_i$ctval
      
      f_kde=NULL
      f_kde=all_f_kde[sample_size_i] #get the kde's sample sizes 
      
      fit <- bbmle::mle2(prev_likelihood_real,start=list(ini_p=runif(1)),
                         data=list(f_kde=f_kde,n_pools=b_pools,n_sample_size=sample_size_i,
                                   observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
                         method="Brent",lower=0,upper=1)
      suppressMessages(
        confint1 <- tryCatch({
          ci <- confint(fit)
          if(is.na(ci[1])) ci[1] <- 0
          if(is.na(ci[2])) ci[1] <- 1
          ci
        }, error = function(e) {
          estms<-sapply(1:length(p_tests),function(x) prev_likelihood_real(f_kde=f_kde,ini_p=p_tests[x],n_pools=b_pools,
                                                                           n_sample_size=sample_size_i,
                                                                           observed_ct_vec=observed_ct_vec,neg_lik=FALSE))
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
          estms<-sapply(1:length(p_tests),function(x) prev_likelihood_real(f_kde=f_kde,ini_p=p_tests[x],n_pools=b_pools,
                                                                           n_sample_size=sample_size_i,
                                                                           observed_ct_vec=observed_ct_vec,neg_lik=FALSE))
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
      
      # optim_prev<-optim(runif(1),fn=prev_likelihood_real,f_kde=f_kde,n_pools=b_pools,n_sample_size=sample_size_i,
      # observed_ct_vec=observed_ct_vec,neg_lik=TRUE,method = "Brent",lower=0,upper=1)$par
      optim_prev_esti<-rbind(optim_prev_esti,data.frame("yer"=yr_to_run,"week"=i,"prev_esti"=coef(fit), "q1"=confint1[1],
                                                        "q2"=confint1[2]))
      
      
    }else{
      optim_prev_esti=data.frame("yer"=yr_to_run,"week"=i,"prev_esti"=NA,"q1"=NA,"q2"=NA)
    }
    
    prev_esti_by_ys_neb=rbind(prev_esti_by_ys_neb,optim_prev_esti)
  }
}

#p_tests<-seq(0,.1,0.002)
#repeat the procedure for each year: for colorado 
prev_esti_by_ys_col<-NULL

for (j in 1:3) {
  yr_to_run=yrs[j]
  cts_of_yr<-subset(dat_col,surv_year==yr_to_run)
  optim_prev_esti<-rep(NA,52)
  #run the prevalence estimation in parallel:
  optim_prev_esti<-NULL
  for(i in 1:52) {
    cts_at_i<-subset(cts_of_yr,cts_of_yr$disease_week==i)
    #remove the row with NA ct values 
    cts_at_i=cts_at_i[complete.cases(cts_at_i),]
    b_pools=nrow(cts_at_i)
    if(b_pools>0){
      sample_size_i=cts_at_i$pool_size
      observed_ct_vec=cts_at_i$ctval
      
      f_kde=NULL
      f_kde=all_f_kde[sample_size_i] #get the kde's sample sizes 
      
      fit <- bbmle::mle2(prev_likelihood_real,start=list(ini_p=runif(1)),
                         data=list(f_kde=f_kde,n_pools=b_pools,n_sample_size=sample_size_i,
                                   observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
                         method="Brent",lower=0,upper=1)
      suppressMessages(
        confint1 <- tryCatch({
          ci <- confint(fit)
          if(is.na(ci[1])) ci[1] <- 0
          if(is.na(ci[2])) ci[1] <- 1
          ci
        }, error = function(e) {
          estms<-sapply(1:length(p_tests),function(x) prev_likelihood_real(f_kde=f_kde,ini_p=p_tests[x],n_pools=b_pools,
                                                                           n_sample_size=sample_size_i,
                                                                           observed_ct_vec=observed_ct_vec,neg_lik=FALSE))
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
          estms<-sapply(1:length(p_tests),function(x) prev_likelihood_real(f_kde=f_kde,ini_p=p_tests[x],n_pools=b_pools,
                                                                           n_sample_size=sample_size_i,
                                                                           observed_ct_vec=observed_ct_vec,neg_lik=FALSE))
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
      
      # optim_prev<-optim(runif(1),fn=prev_likelihood_real,f_kde=f_kde,n_pools=b_pools,n_sample_size=sample_size_i,
      # observed_ct_vec=observed_ct_vec,neg_lik=TRUE,method = "Brent",lower=0,upper=1)$par
      optim_prev_esti<-rbind(optim_prev_esti,data.frame("yer"=yr_to_run,"week"=i,"prev_esti"=coef(fit), "q1"=confint1[1],
                                                        "q2"=confint1[2]))
      
    }else{
      optim_prev_esti=data.frame("yer"=yr_to_run,"week"=i,"prev_esti"=NA,"q1"=NA,"q2"=NA)
    }
    prev_esti_by_ys_col=rbind(prev_esti_by_ys_col,optim_prev_esti)
  }
}

prev_esti_by_ys_neb$state<-rep("Nebraska")
prev_esti_by_ys_col$state<-rep("Colorado")

all_esti<-rbind(prev_esti_by_ys_neb,prev_esti_by_ys_col)


suum<-all_esti %>%
  group_by(state,yer) %>%
  summarise(peak=max(prev_esti,na.rm = T),max_q1=q1[which.max(prev_esti)],max_q2=q2[which.max(prev_esti)])
suum
#using the cdc package 


pv_cdc_neb<-data.frame(matrix(NA,nrow = 0,ncol = 6))
colnames(pv_cdc_neb)<-c("state","year","week","cdc_prev","q1","q2")

#calculate the prevalence from the cdc package
for (j in 1:3) {
  yr_to_run=yrs[j]
  cts_of_yr<-subset(dat_neb,surv_year==yr_to_run)
  for(i in 1:52){
    cts_at_i<-subset(cts_of_yr,cts_of_yr$disease_week==i)
    #remove the row with NA ct values 
    cts_at_i=cts_at_i[complete.cases(cts_at_i),]
    #inf_states<-(cts_at_i$ctval<35)*1
    inf_states<-cts_at_i$WNV
    pol_sz<-cts_at_i$pool_size
    if(sum(pol_sz)>0){
      pv<-pooledBin(inf_states,pol_sz,ci.method ="lrt" )
      pv_cdc_neb<-rbind(pv_cdc_neb,data.frame("state"="Nebraska","year"=yr_to_run,"week"=i,"cdc_prev"=pv$P,"q1"=pv$Lower,"q2"=pv$Upper))
    }
    else{
      pv_cdc_neb<-rbind(pv_cdc_neb,data.frame("state"="Nebraska","year"=yr_to_run,"week"=i,"cdc_prev"=NA,"q1"=NA,"q2"=NA))
    }
  }
}

pv_cdc_col<-data.frame(matrix(NA,nrow = 0,ncol = 6))
colnames(pv_cdc_col)<-c("state","year","week","cdc_prev","q1","q2")

#calculate the prevalence from the cdc package
for (j in 1:3) {
  yr_to_run=yrs[j]
  cts_of_yr<-subset(dat_col,surv_year==yr_to_run)
  for(i in 1:52){
    cts_at_i<-subset(cts_of_yr,cts_of_yr$disease_week==i)
    #remove the row with NA ct values 
    cts_at_i=cts_at_i[complete.cases(cts_at_i),]
    #inf_states<-(cts_at_i$ctval<40)*1
    inf_states<-cts_at_i$WNV
    pol_sz<-cts_at_i$pool_size
    if(sum(pol_sz)>0){
      pv<-pooledBin(inf_states,pol_sz,ci.method ="lrt" )
      pv_cdc_col<-rbind(pv_cdc_col,data.frame("state"="Colorado","year"=yr_to_run,"week"=i,"cdc_prev"=pv$P,"q1"=pv$Lower,"q2"=pv$Upper))
    }
    else{
      pv_cdc_col<-rbind(pv_cdc_col,data.frame("state"="Colorado","year"=yr_to_run,"week"=i,"cdc_prev"=NA,"q1"=NA,"q2"=NA))
    }
  }
}

all_cdc<-rbind(pv_cdc_neb,pv_cdc_col)

suum<-all_cdc %>%
  group_by(state,yer) %>%
  summarise(peak=max(cdc_prev,na.rm = T),max_q1=q1[which.max(cdc_prev)],max_q2=q2[which.max(cdc_prev)])
suum



coln<-c("yer", "week" , "prev","q1", "q2","state","method")
all_esti$method=rep("pooled",nrow(all_esti))
colnames(all_esti)=coln
all_cdc$method=rep("cdc")
a_cdc=data.frame(all_cdc$year,all_cdc$week,all_cdc$cdc_prev,all_cdc$q1,all_cdc$q2,all_cdc$state,all_cdc$method)
colnames(a_cdc)=coln

all_prev_est<-rbind(all_esti,a_cdc)



p1<-ggplot(data=subset(all_prev_est,state=="Colorado"),aes(x=week,y=prev))+
  geom_ribbon(aes(x=week,ymin=q1,ymax=q2,fill = method),alpha=1)+
  geom_line()+
  facet_wrap(as.factor(yer)~method,ncol=6)+
  scale_fill_manual(values = c("#a86463","#ad7f00"),
                    labels=c("Binarised data","Pooled Ct values"))+
  xlim(20,45)+
  ylim(0,0.025)+
  ylab("Estimated Prevalence")+
  xlab("Week")+
  scale_color_manual(values  = wes_palette("AsteroidCity1", 3, type = "continuous"))+
  default_theme +
  ggtitle("(B) Colorado")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p1

ggsave("figures/p1.png",p1,height = 8,width = 14)



p2<-ggplot(data=subset(all_prev_est,state=="Nebraska"),aes(x=week,y=prev))+
  geom_ribbon(aes(x=week,ymin=q1,ymax=q2,fill = method),alpha=1)+
  geom_line()+
  facet_wrap(~as.factor(yer)~method,ncol=6)+
  scale_fill_manual(values = c("#a86463","#ad7f00"),
                    labels=c("Binarised data","Pooled Ct values"))+
  xlim(20,45)+
  # ylim(0,0.08)+
  ylab("Estimated Prevalence")+
  xlab("Week")+
  scale_color_manual(values  = wes_palette("AsteroidCity1", 3, type = "continuous"))+
  default_theme +
  ggtitle("(A) Nebraska")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p2


ggarrange(p2,p1,nrow = 2,common.legend = T, legend="bottom")

#ggsave("figures/figure_7.png",last_plot(),height = 8,width = 14)
