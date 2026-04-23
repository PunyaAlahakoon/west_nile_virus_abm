##estimating the productive vs non-productive prevalance with data


#load the default theme for figures:
source("functions/default_theme.R")
default_theme<-default_theme()

#read all the fkdes 
all_f_kde=NULL
all_f_kde=sapply(1:50,function(x) readRDS(file = paste0("pre_calculations/kde_calcs/f_kde_",x,".rds")))


#joint KDEs taht were precalculated. See KDE calculations to see how they were pre-caluculated
all_f_kde_joint=NULL
all_f_kde_joint=lapply(1:50,function(x) readRDS(file = paste0("pre_calculations/joint_kde_calc/joint_fkde_",x,".rds")))


#read the clean real datasets 
source("functions/clean_nebraska_data.R")
dat_neb<-clean_nebraska_data()
source("functions/clean_colorado_data.R")
dat_col<-clean_colorado_data()

#likelihood functions 
source("functions/prev_likelihood_real_static_dynamic5_vec_reparametrise2.R")
source("functions/prev_likelihood_real.R")



###additional estimates 
p_tests<-seq(0,.5,0.005)

#extract the disease weeks and years 
yrs<-2022:2024
weeks<-24:39


#precompute the possible combinations mat; so that joint KDE combinations are the same as this: 
# Initialize list to store combinations for each pool size
all_combos <- vector("list", 50)
for (n in 1:50) {
  # For each pool size n, generate all combinations of (static, dynamic, negatives)
  combos_n <- expand.grid(
    static = 0:n,
    dynamic = 0:n,
    negatives = 0:n
  )
  
  # Keep only rows where sum == n
  combos_n <- subset(combos_n, static + dynamic + negatives == n)
  
  # Store in list
  all_combos[[n]] <- combos_n
}



S <- 5000 #for bootstrapping 

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
          if(is.na(ci[2])) ci[2] <- 1
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
      
      
      f_kde_joint=NULL
      #  f_kde_joint=all_f_kde_joint[sample_size_i] #get the kde's sample sizes 
      
      f_kde_joint=all_f_kde_joint #get the kde's sample sizes 
      
      fit_a <- bbmle::mle2(prev_likelihood_real_static_dynamic5_vec_reparametrise2,
                           start=list(a = qlogis(0.005), b = qlogis(0.005)),
                           data=list(f_kde_joint=f_kde_joint,all_combos=all_combos,
                                     n_pools=b_pools,
                                     pool_sizes=sample_size_i,
                                     observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
                           method="L-BFGS-B",
                           lower=c(-50,-50),upper=c(50,50))
      
      a_hat <- coef(fit_a)[["a"]]
      b_hat <- coef(fit_a)[["b"]]
      p_static_hat <- plogis(a_hat)
      p_dynamic_hat <- plogis(b_hat) * (1 - p_static_hat)
      
      boots <- MASS::mvrnorm(S, mu = coef(fit_a), Sigma = vcov(fit_a))
      a_sim <- boots[,1]
      b_sim <- boots[,2]
      
      p_static_sim  <- plogis(a_sim)
      p_dynamic_sim <- plogis(b_sim) * (1 - p_static_sim)
      #also calculate the proportions:
      p_dynamic_sim_prop<-p_dynamic_sim/(p_dynamic_sim+p_static_sim)
      p_static_sim_prop<-p_static_sim/(p_dynamic_sim+p_static_sim)
      
      #also calculate the total prevalence<
      p_calc_total<-p_dynamic_sim+p_static_sim
      
      #confidence intervals 
      ci_p_static  <- quantile(p_static_sim,  probs = c(.025, .975))
      ci_p_dynamic <- quantile(p_dynamic_sim, probs = c(.025, .975))
      
      ci_prop_static<-quantile(p_static_sim_prop,  probs = c(.025, .975),na.rm=T)
      ci_prop_dynamic<-quantile(p_dynamic_sim_prop,  probs = c(.025, .975),na.rm=T)
      
      ci_total_calc<-quantile(p_calc_total,probs = c(.025, .975),na.rm=T)
      
      # optim_prev<-optim(runif(1),fn=prev_likelihood_real,f_kde=f_kde,n_pools=b_pools,n_sample_size=sample_size_i,
      # observed_ct_vec=observed_ct_vec,neg_lik=TRUE,method = "Brent",lower=0,upper=1)$par
      optim_prev_esti<-rbind(optim_prev_esti,data.frame(yer=yr_to_run,
                                                        week=i,
                                                        prev_esti=coef(fit), 
                                                        q1=confint1[1],
                                                        q2=confint1[2],
                              optim_p_static   =p_static_hat,
                             optim_p_dynamic=p_dynamic_hat,
                             static_q1      = ci_p_static[1],
                             static_q2      = ci_p_static[2],
                             dynamic_q1      = ci_p_dynamic[1],
                             dynamic_q2      = ci_p_dynamic[2],
                             prop_dynamic_q1=ci_prop_dynamic[1],
                             prop_dynamic_q2=ci_prop_dynamic[2],
                             prop_static_q1=ci_prop_static[1],
                             prop_static_q2=ci_prop_static[2],
                             ci_total_calc_q1=ci_total_calc[1],
                             ci_total_calc_q2=ci_total_calc[2]
      )
                             )
      
      
    }else{
      optim_prev_esti=data.frame(yer=yr_to_run,
                                 week=i,
                                 prev_esti=NA,
                                 q1=NA,
                                 q2=NA,
                                 optim_p_static=NA,
                                 optim_p_dynamic=NA,
                                 static_q1      = NA,
                                 static_q2      = NA,
                                 dynamic_q1      = NA,
                                 dynamic_q2      = NA,
                                 prop_dynamic_q1=NA,
                                 prop_dynamic_q2=NA,
                                 prop_static_q1=NA,
                                 prop_static_q2=NA,
                                 ci_total_calc_q1=NA,
                                 ci_total_calc_q2=NA
                                 )
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
          if(is.na(ci[2])) ci[2] <- 1
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
      
      
      
      f_kde_joint=NULL
      #  f_kde_joint=all_f_kde_joint[sample_size_i] #get the kde's sample sizes 
      
      f_kde_joint=all_f_kde_joint #get the kde's sample sizes 
      
      fit_a <- bbmle::mle2(prev_likelihood_real_static_dynamic5_vec_reparametrise2,
                           start=list(a = qlogis(0.005), b = qlogis(0.005)),
                           data=list(f_kde_joint=f_kde_joint,all_combos=all_combos,
                                     n_pools=b_pools,
                                     pool_sizes=sample_size_i,
                                     observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
                           method="L-BFGS-B",
                           lower=c(-50,-50),upper=c(50,50))
      
      a_hat <- coef(fit_a)[["a"]]
      b_hat <- coef(fit_a)[["b"]]
      p_static_hat <- plogis(a_hat)
      p_dynamic_hat <- plogis(b_hat) * (1 - p_static_hat)
      
      boots <- MASS::mvrnorm(S, mu = coef(fit_a), Sigma = vcov(fit_a))
      a_sim <- boots[,1]
      b_sim <- boots[,2]
      
      p_static_sim  <- plogis(a_sim)
      p_dynamic_sim <- plogis(b_sim) * (1 - p_static_sim)
      #also calculate the proportions:
      p_dynamic_sim_prop<-p_dynamic_sim/(p_dynamic_sim+p_static_sim)
      p_static_sim_prop<-p_static_sim/(p_dynamic_sim+p_static_sim)
      
      #also calculate the total prevalence<
      p_calc_total<-p_dynamic_sim+p_static_sim
      
      #confidence intervals 
      ci_p_static  <- quantile(p_static_sim,  probs = c(.025, .975))
      ci_p_dynamic <- quantile(p_dynamic_sim, probs = c(.025, .975))
      
      ci_prop_static<-quantile(p_static_sim_prop,  probs = c(.025, .975),na.rm=T)
      ci_prop_dynamic<-quantile(p_dynamic_sim_prop,  probs = c(.025, .975),na.rm=T)
      
      ci_total_calc<-quantile(p_calc_total,probs = c(.025, .975),na.rm=T)
      
      
      
      
      # optim_prev<-optim(runif(1),fn=prev_likelihood_real,f_kde=f_kde,n_pools=b_pools,n_sample_size=sample_size_i,
      # observed_ct_vec=observed_ct_vec,neg_lik=TRUE,method = "Brent",lower=0,upper=1)$par
      optim_prev_esti<-rbind(optim_prev_esti,data.frame("yer"=yr_to_run,
                                                        "week"=i,
                                                        "prev_esti"=coef(fit),
                                                        "q1"=confint1[1],
                                                        "q2"=confint1[2],
                                                        optim_p_static   =p_static_hat,
                                                        optim_p_dynamic=p_dynamic_hat,
                                                        static_q1      = ci_p_static[1],
                                                        static_q2      = ci_p_static[2],
                                                        dynamic_q1      = ci_p_dynamic[1],
                                                        dynamic_q2      = ci_p_dynamic[2],
                                                        prop_dynamic_q1=ci_prop_dynamic[1],
                                                        prop_dynamic_q2=ci_prop_dynamic[2],
                                                        prop_static_q1=ci_prop_static[1],
                                                        prop_static_q2=ci_prop_static[2],
                                                        ci_total_calc_q1=ci_total_calc[1],
                                                        ci_total_calc_q2=ci_total_calc[2]
                                                        )
                             
                             )
      
    }else{
      optim_prev_esti=data.frame("yer"=yr_to_run,
                                 "week"=i,
                                 "prev_esti"=NA,
                                 "q1"=NA,
                                 "q2"=NA,
                                 "optim_p_static"=NA,
                                 "optim_p_dynamic"=NA,
                                 static_q1      = NA,
                                 static_q2      = NA,
                                 dynamic_q1      = NA,
                                 dynamic_q2      = NA,
                                 prop_dynamic_q1=NA,
                                 prop_dynamic_q2=NA,
                                 prop_static_q1=NA,
                                 prop_static_q2=NA,
                                 ci_total_calc_q1=NA,
                                 ci_total_calc_q2=NA
                                 )
    }
    prev_esti_by_ys_col=rbind(prev_esti_by_ys_col,optim_prev_esti)
  }
}

prev_esti_by_ys_neb$state<-rep("Nebraska")
prev_esti_by_ys_col$state<-rep("Colorado")



pv_cdc_neb<-data.frame(matrix(NA,nrow = 0,ncol = 6))
colnames(pv_cdc_neb)<-c("state","yer","week","cdc_prev","q1","q2")

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
      pv_cdc_neb<-rbind(pv_cdc_neb,data.frame("state"="Nebraska","yer"=yr_to_run,"week"=i,"cdc_prev"=pv$P,"q1"=pv$Lower,"q2"=pv$Upper))
    }
    else{
      pv_cdc_neb<-rbind(pv_cdc_neb,data.frame("state"="Nebraska","yer"=yr_to_run,"week"=i,"cdc_prev"=NA,"q1"=NA,"q2"=NA))
    }
  }
}

pv_cdc_col<-data.frame(matrix(NA,nrow = 0,ncol = 6))
colnames(pv_cdc_col)<-c("state","yer","week","cdc_prev","q1","q2")

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
      pv_cdc_col<-rbind(pv_cdc_col,data.frame("state"="Colorado","yer"=yr_to_run,"week"=i,"cdc_prev"=pv$P,"q1"=pv$Lower,"q2"=pv$Upper))
    }
    else{
      pv_cdc_col<-rbind(pv_cdc_col,data.frame("state"="Colorado","yer"=yr_to_run,"week"=i,"cdc_prev"=NA,"q1"=NA,"q2"=NA))
    }
  }
}

all_cdc<-rbind(pv_cdc_neb,pv_cdc_col)


all_esti<-rbind(prev_esti_by_ys_neb,prev_esti_by_ys_col)

all_esti_unique <- unique(all_esti)
#all_esti_unique=all_esti_unique[-1611,]


library(dplyr)


# str(all_cdc[, c("yer", "week", "state")])
# str(all_esti_unique[, c("yer", "week", "state")])


combined <- all_cdc %>%
  left_join(all_esti_unique,
            by = c("yer" = "yer", "week" = "week", "state" = "state"),
            suffix = c("_cdc", "_esti"))

#saveRDS(combined,"combined.RDS")

combined <- readRDS("combined.RDS")
pl_all_esti<-combined[,c("yer","week","cdc_prev","prev_esti","optim_p_dynamic","state")]
#pl_all_esti<-combined[,c(2:3, 4,7, 11,1)]
colnames(pl_all_esti)<-c("yer","week", "cdc", "overall", "dynamic",'state')
m_pl_all_esti<-melt(pl_all_esti,id=c("yer","week",'state'))

all_esti_q1s<-combined[,c("yer","week","q1_cdc","q1_esti","dynamic_q1","state")]
colnames(all_esti_q1s)<-c("yer","week","cdc", "overall", "dynamic",'state')
all_esti_q2s<-combined[,c("yer","week","q2_cdc","q2_esti","dynamic_q2","state")]
colnames(all_esti_q2s)<-c("yer","week","cdc", "overall", "dynamic",'state')

m_all_esti_q1s<-melt(all_esti_q1s,id=c("yer","week",'state'))
names(m_all_esti_q1s)[names(m_all_esti_q1s) == "value"] <- "q1"

m_all_esti_q2s<-melt(all_esti_q2s,id=c("yer","week",'state'))
names(m_all_esti_q2s)[names(m_all_esti_q2s) == "value"] <- "q2"
m_all_esti_qs<-cbind(m_all_esti_q1s,"q2"=m_all_esti_q2s$q2)


est_true_prev<-subset(m_pl_all_esti,state=="Nebraska" & !variable %in% c("overall","dynamic")) #just keep cdc
levels(est_true_prev$variable)<-droplevels(est_true_prev$variable)
levels(est_true_prev$variable)<-"dynamic" #to plot this over the dynamic 
est_true_prev_2<-est_true_prev
est_true_prev_2$variable<-rep("overall")
est_true_prev<-rbind(est_true_prev,est_true_prev_2)
  
p1<-ggplot(data=subset(m_pl_all_esti,state=="Nebraska"))+
  geom_ribbon(data=subset(m_all_esti_qs,state=="Nebraska"),
              aes(x=week,ymin = q1,ymax = q2,fill =variable ),alpha=0.6)+
  geom_line(aes(x=week,y=value,color=variable))+
  scale_fill_manual(values = c("#a86463","#5d861f","#ad7f00"),
                    labels=c("Binarised data-based","Pooled Ct value-based ",
                             "Pooled Ct value-based \n (Productive infection)"))+
  scale_color_manual(values = c("#971a18","#385905","#674D06") )+
  geom_line(data=est_true_prev,aes(x=week,y=value),color="#971a18",linetype="dashed")+
  coord_cartesian(ylim=c(0,0.04))+
  xlim(20,40)+
  ylab(NULL)+
  # 
  facet_grid(
    rows = vars(yer),
    cols = vars(variable),
    switch = "y"
  ) +
  
  default_theme+
  ggtitle("(A) Nebraska")+
  theme(
    strip.background = element_blank(),
    
    # Remove variable labels
    strip.text.x = element_blank(),
    

    strip.placement = "outside",
    strip.text.y.right = element_text(size = 12),
    
    legend.position = "bottom"
  ) +
  guides(color = "none")

p1

est_true_prev<-subset(m_pl_all_esti,state=="Colorado" & !variable %in% c("overall","dynamic") & yer!=2024 )
levels(est_true_prev$variable)<-droplevels(est_true_prev$variable)
levels(est_true_prev$variable)<-"dynamic"
est_true_prev_2<-est_true_prev
est_true_prev_2$variable<-rep("overall")
est_true_prev<-rbind(est_true_prev,est_true_prev_2)

p2 <- ggplot(data=subset(m_pl_all_esti, state=="Colorado" & yer != 2024)) +
  geom_ribbon(
    data = subset(m_all_esti_qs, state=="Colorado" & yer != 2024),
    aes(x = week, ymin = q1, ymax = q2, fill = variable),
    alpha = 0.7
  ) +
  geom_line(aes(x=week,y=value,color=variable))+
  scale_fill_manual(values = c("#a86463","#5d861f","#ad7f00"),
                    labels=c("Binarised data-based","Pooled Ct value-based ",
                             "Pooled Ct value-based \n (Productive infection)"))+
  scale_color_manual(values = c("#971a18","#385905","#674D06") )+
  geom_line(data=est_true_prev,aes(x=week,y=value),color="#971a18",linetype="dashed")+
  geom_line(
    data = est_true_prev,
    aes(x = week, y = value),
    color = "#971a18",
    linetype = "dashed"
  ) +
  coord_cartesian(ylim = c(0,0.02)) +
  xlim(20,40) +
  ylab(NULL)+
  # 
  facet_grid(
    rows = vars(yer),
    cols = vars(variable),
    switch = "y"
  ) +
  
  default_theme +
  ggtitle("(B) Colorado") +
  theme(
    strip.background = element_blank(),
  
    strip.text.x = element_blank(),
    
    strip.placement = "outside",
    strip.text.y.right = element_text(size = 12),
    
    legend.position = "bottom"
  ) +
  guides(color = "none")

p2

fig <- ggarrange(
  p1 + coord_cartesian(ylim=c(0, 0.04)),
  p2 + coord_cartesian(ylim = c(0, 0.02)),
  nrow = 2,
  align = "v",
  common.legend = T,
  legend = "bottom",
  heights = c(2, 1.5)
)

annotate_figure(fig, left = text_grob("Estimated prevalence", rot = 90,size=14))

#ggarrange(p1,p2,nrow = 2)

ggsave("figures/figure_7_v2.png",last_plot(),height = 11,width = 12)


# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(forcats)
# 


combined$prop_dynamic<-combined$optim_p_dynamic/(combined$optim_p_dynamic+combined$optim_p_static)
combined$prop_static<-combined$optim_p_static/(combined$optim_p_dynamic+combined$optim_p_static)

df <- combined %>%
  mutate(
    week = as.numeric(week),
    yer  = as.factor(yer),
    state = as.factor(state)
  )


long_ribbons <- df %>%
  transmute(
    state, yer, week,
    system = "Dynamic",
    q1 = prop_dynamic_q1,
    q2 = prop_dynamic_q2,
     est = prop_dynamic,
   
   )

all_dat<-rbind(dat_neb,dat_col)
names(all_dat)[names(all_dat) == "surv_year"] <- "yer"
names(all_dat)[names(all_dat) == "disease_week"] <- "week"

long_ribbons_aug <- long_ribbons %>%
  left_join(
    all_dat %>% select(state, yer, week, ctval, pool_size),
    by = c("state", "yer", "week")
  )
  


#%>%
  # bind_rows(
  #   df %>%
  #     transmute(
  #       state, yer, week,
  #       system = "Static",
  #       q1 = prop_static_q1,
  #       q2 = prop_static_q2,
  # 
  #       est = optim_p_static/(optim_p_dynamic+optim_p_static),
  #     )
  # ) %>%
  # mutate(system = factor(system, levels = c("Static", "Dynamic")))  # set order in legend


long_ribbons <- long_ribbons %>%
  filter(!(is.na(q1) & is.na(q2)))

long_ribbons$state=factor(long_ribbons$state)
long_ribbons$state=factor(long_ribbons$state,levels=c("Nebraska","Colorado"))

p <- ggplot(long_ribbons, aes(x = week, group = yer)) +
  # Ribbon for the interval
  geom_ribbon(
    aes(ymin = q1, ymax = q2, fill = system),
    alpha = 0.20,
    colour = NA
  ) +
  geom_line(
    aes(y = est, colour = system),
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  facet_wrap(state~yer, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("Static" = "#1f77b4", "Dynamic" = "darkred")) +
  scale_fill_manual(values   = c("Static" = "#1f77b4", "Dynamic" = "darkred")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_hline(yintercept = 0.5, color = "black", linetype="dashed")+
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Week",
    y = "Proportion",
    colour = "System",
    fill   = "",
    linetype = "Year",
    title = "Estimated productive infection proportions by state"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

p

ggsave("figures/suplement_prod_prop.png",p,height =6,width = 12)

# 
# long_ribbons_aug <- long_ribbons_aug %>%
#   filter(!(is.na(ctval)))
# 
# 
# counts_yw <- long_ribbons_aug %>%
#   count(yer, week,state)
# head(counts_yw)
# 
# 
# 
# long_ribbons_x <- long_ribbons %>%
#   left_join(
#     counts_yw %>% select(state, yer, week, n),
#     by = c("state", "yer", "week")
#   )
# 
# head(long_ribbons_x)
# 
# p_bottom <- ggplot(counts_yw, aes(x = week, y= n)) +
#   geom_col(fill = "grey50", width = 0.8) +
#   facet_wrap(yer ~ state, scales = "free_y", ncol = 3) +
#   labs(y = "Pools", x = "Week") +
#   theme_minimal()
# 
# p_bottom
# 
# p_cts<-ggplot(long_ribbons_aug, aes(x = week, y = ctval)) +
#   geom_point()+
#   facet_wrap(yer ~ state, scales = "free_y", ncol = 3) +
#   labs(y = "Ct", x = "Week") +
#   theme_minimal()
# 
# 
# pool_size<-ggplot(long_ribbons_aug, aes(x = week, y =pool_size )) +
#   geom_point(aes(colour =(pool_size)))+
#   facet_wrap(yer ~ state, scales = "free_y", ncol = 3) +
#   labs(y = "Ct", x = "Week") +
#   theme_minimal()
# pool_size
# 
# p/p_cts/pool_size
# 
# 
# ggplot(data=long_ribbons_aug,aes(x=week,y=pool_size,fill=ctval))+
#   geom_tile()+  
#   facet_grid(yer~state)
# 
# #number of pools that are positive 
# counts_yw <- subset(long_ribbons_aug,ctval<40) %>%
#   count(yer, week,state)
# head(counts_yw)
# 
# #proportion:
# 
# counts_yw <- long_ribbons_aug %>%
#   group_by(state, week,yer) %>%
#   summarise(
#     prop_ct_lt40 = mean(ctval < 40, na.rm = TRUE)
#   )
# counts_yw
# 
# ggplot(counts_yw, aes(x = week, y= prop_ct_lt40)) +
#   geom_line() +
#   facet_wrap(yer ~ state, scales = "free_y", ncol = 3) +
#   labs(y = "Positive pools", x = "Week") +
#   theme_minimal()
# 
# 
# px <- ggplot(long_ribbons_x, aes(x = week, y =est)) +
# geom_line()+
#   facet_wrap(yer~ state, scales = "free_y", ncol = 3) 
#  
# px
# 
# 
# 
# 
# coeff<-10
# p <- ggplot(long_ribbons, aes(x = week, group = yer)) +
#   # Ribbon for the interval
#   geom_ribbon(
#     aes(ymin = q1, ymax = q2, fill = system),
#     alpha = 0.20,
#     colour = NA
#   ) +
#   geom_line(
#     aes(y = est, colour = system),
#     linewidth = 0.9,
#     na.rm = TRUE
#   ) +
#   geom_line(data=counts_yw,aes(x=week,y=-prop_ct_lt40*coeff))+
#   scale_y_continuous(
#     
#     # Features of the first axis
#     name = "Proportion",
#     
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~-./coeff, name="Productive infection prevalance")
#   )+
#   facet_wrap(yer~ state, scales = "free_y", ncol = 3) +
#   scale_colour_manual(values = c("Static" = "#1f77b4", "Dynamic" = "#d62728")) +
#   scale_fill_manual(values   = c("Static" = "#1f77b4", "Dynamic" = "#d62728")) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
#   geom_hline(yintercept = 0.5, color = "blue")+
#   # #scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   # labs(
#   #   x = "Week",
#   #   y = "Proportion",
#   #   colour = "System",
#   #   fill   = "",
#   #   linetype = "Year",
#   #   title = "Proportions by State",
#   #   subtitle = "Static vs Dynamic"
#   # ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     legend.position = "none",
#     panel.grid.minor = element_blank()
#   )
# 
# p
# 
# 
# avg_by_week <- long_ribbons %>%
#   group_by(week,state) %>%
#   summarise(
#     mean_est = mean(est, na.rm = TRUE),
#     n = sum(!is.na(est))
#   ) %>%
#   arrange(week)
# 
# avg_by_week
# 
# 
# 
# ggplot(data=avg_by_week)+
#   geom_line(aes(x=week,y=mean_est))+
#   geom_hline(yintercept = 0.5, color = "blue")+
#   facet_grid(~state)
# 
# 
# 
# 
# x1<-ggplot(combined)+
#   geom_line(aes(x=week,y=prev_esti),color="black")+
#   geom_line(aes(x=week,y=optim_p_static+optim_p_dynamic),color="red")+
#   facet_wrap(state~yer)
# x1
# 
# x2<-ggplot(combined)+
#   geom_line(aes(x=week,y=(optim_p_dynamic/(optim_p_dynamic+optim_p_static))))+
#   facet_wrap(state~yer)
# x2
# 
# x1/x2
# 
# 
# 
# 
# 
# 
# peaks <- combined %>%
#   group_by(state, yer) %>%
#   slice_max(order_by = optim_p_dynamic, n = 1, with_ties = FALSE) 
# peaks
# 
# 
# p1<-ggplot(data=subset(m_pl_all_esti,state=="Nebraska"))+
#   geom_ribbon(data=subset(m_all_esti_qs,state=="Nebraska"),
#               aes(x=week,ymin = q1,ymax = q2,fill =variable ))+
#   geom_line(aes(x=week,y=value))+
#   scale_fill_manual(values = c("#a86463","#ad7f00","#5d861f"),
#                      labels=c("Pooled Ct value-based ","Binarised data-based",
#                               "Pooled Ct value-based \n (Productive infection)"))+
#                 
#   ylim(0,0.05)+
#   facet_wrap(variable~yer)+
#   default_theme+
#   ggtitle("(A) Nebraska")+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )+
#   theme(legend.position = "none")
# 
# 
# p1
# 
# p2<-ggplot(data=subset(m_pl_all_esti,state=="Colorado" & yer!=2024 ))+
#   geom_ribbon(data=subset(m_all_esti_qs,state=="Colorado" & yer!=2024 ),
#               aes(x=week,ymin = q1,ymax = q2,fill =variable ))+
#   geom_line(aes(x=week,y=value))+
#   scale_fill_manual(values = c("#a86463","#ad7f00","#5d861f"),
#                     labels=c("Pooled Ct value-based ","Binarised data-based",
#                              "Pooled Ct value-based \n (Productive infection)"))+
#   
#   ylim(0,0.05)+
#   facet_wrap(variable~yer,ncol=2)+
#   default_theme+
#   ggtitle("(B) Colorado")+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )
# 
# 
# p2
# 
# 
# 
# ggarrange(p1,p2,nrow = 2)
# 
# ggsave("figures/figure_7_v2.png",last_plot(),height = 14,width = 14)
# 
# 

