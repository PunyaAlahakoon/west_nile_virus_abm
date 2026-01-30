#load the default theme for figures:
theme<-source("functions/default_theme.R")
default_theme<-default_theme()

source("functions/prev_likelihood_real_static_dynamic.R")
source("functions/prev_likelihood_real.R")
#load all the data:
all_data <- readRDS("synthetic_data/all_data.rds")
#load the data: pooled ct values 
cts_per_pool= readRDS("pre_calculations/pooled_cts/cts_per_pool_demo.rds") #they are only data from 2:4 years for the 4 sims, 

ct_threshold=40


#also doubke check if the sample size is greater than 50:
cts_per_pool$sample_size[cts_per_pool$sample_size>50]=50


#KDE estimation:
#positive_ctss<-all_data %>% filter(ct_value < ct_threshold) %>% pull(ct_value)

#choose the sim number and the year to run:
yr_to_run=4
sim_to_run=2
cts_per_poolss<-subset(cts_per_pool,(cts_per_pool$sim_number==sim_to_run) & (cts_per_pool$year==yr_to_run))



#read all the fkdes 
#read all the fkdes 
all_f_kde=NULL
all_f_kde=sapply(1:50,function(x) readRDS(file = paste0("pre_calculations/kde_calcs/f_kde_",x,".rds")))
# 
# all_f_kde_static=NULL
# all_f_kde_static=sapply(1:50,function(x) readRDS(file = paste0("pre_calculations/kde_calcs_dynamic_static/f_kde_static_",x,".rds")))
# 
# all_f_kde_dynamic=NULL
# all_f_kde_dynamic=sapply(1:50,function(x) readRDS(file = paste0("pre_calculations/kde_calcs_dynamic_static/f_kde_dynamic_",x,".rds")))


#calculate the prevalances: 
all_f_kde_joint=NULL
all_f_kde_joint=lapply(1:50,function(x) readRDS(file = paste0("pre_calculations/joint_kde_calc/joint_fkde_",x,".rds")))

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


#for bootstrapping 
set.seed(1)
S <- 2000

#run the prevalence estimation: 
optim_prev_esti=NULL


for (i in 1:52) {
  cts_at_i<-subset(cts_per_poolss,cts_per_poolss$week==i)
  b_pools=nrow(cts_at_i)
  if(b_pools>0){
    sample_size_i=cts_at_i$sample_size
    observed_ct_vec=cts_at_i$pooled_ct_per_pool
    
    
    f_kde=NULL
    f_kde=all_f_kde[sample_size_i] #get the kde's sample sizes 
    
    fit <- bbmle::mle2(prev_likelihood_real,start=list(ini_p=runif(1)),
                       data=list(f_kde=f_kde,n_pools=b_pools,n_sample_size=sample_size_i,
                                 observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
                       method="Brent",lower=0,upper=1)
    suppressMessages(
      confint1 <- tryCatch({
        ci <- bbmle::confint(fit)
        if(is.na(ci[1])) ci[1] <- 0
        if(is.na(ci[2])) ci[1] <- 1
        ci
      }, error = function(e) {
        c(0, 1)
      }, warning = function(w) {
        c(0, 1)
      })
    )
    
    #to calculate the the dynamc vs static proposrions, do it independently: 
    # fit_a <- bbmle::mle2(prev_likelihood_real_static_dynamic5_vec_reparametrise2,
    #                      start=list(a = qlogis(0.005), b = qlogis(0.005)),
    #                      data=list(f_kde_joint=all_f_kde_joint,all_combos=all_combos,
    #                                n_pools=b_pools,
    #                                pool_sizes=sample_size_i,
    #                                observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
    #                      method="L-BFGS-B",
    #                      lower=c(-50,-50),upper=c(50,50))
    
    fit_a <- bbmle::mle2(prev_likelihood_real_static_dynamic5_vec_reparametrise2,
                         start=list(a = runif(1), b = runif(1)),
                         data=list(f_kde_joint=all_f_kde_joint,all_combos=all_combos,
                                   n_pools=b_pools,
                                   pool_sizes=sample_size_i,
                                   observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
                         method="L-BFGS-B",
                         lower=c(0,0),upper=c(1,1))
    
    a_hat <- coef(fit_a)[["a"]]
    b_hat <- coef(fit_a)[["b"]]
    p_static_hat <- plogis(a_hat)
    p_dynamic_hat <- plogis(b_hat) * (1 - p_static_hat)
    # 

   
    boots <- MASS::mvrnorm(S, mu = coef(fit_a), Sigma = vcov(fit_a))
    a_sim <- boots[,1]
    b_sim <- boots[,2]
    
    p_static_sim  <- plogis(a_sim)
    p_dynamic_sim <- plogis(b_sim) * (1 - p_static_sim)
    
    ci_p_static  <- quantile(p_static_sim,  probs = c(.025, .975))
    ci_p_dynamic <- quantile(p_dynamic_sim, probs = c(.025, .975))
    
    
    # suppressMessages({
    #   confint_list <- tryCatch({
    #     # Profile likelihood CIs for both parameters
    #     ci_all <- bbmle::confint(fit_a, method = "quad")
    #     #or use this for a a more accurate but reduced number of steps
    #     #using the profilinh method:
    #    # ci_all <- bbmle::confint(fit_a, method = "profile", nsteps = 50)
    #     # Extract as numeric vectors
    #     ci_p_static <- as.numeric(ci_all["p_static", ])
    #     ci_p_dynamic <- as.numeric(ci_all["p_dynamic", ])
    # 
    #     if (is.na(ci_p_static[1])) ci_p_static[1] <- 0
    #     if (is.na(ci_p_static[2])) ci_p_static[2] <- 1
    #     if (is.na(ci_p_dynamic[1])) ci_p_dynamic[1] <- 0
    #     if (is.na(ci_p_dynamic[2])) ci_p_dynamic[2] <- 1
    # 
    #     # Ensure lower <= upper
    #     if (ci_p_static[1] > ci_p_static[2]) ci_p_static <- sort(ci_p_static)
    #     if (ci_p_dynamic[1] > ci_p_dynamic[2]) ci_p_dynamic <- sort(ci_p_dynamic)
    # 
    #     # Return as a named list
    #     list(p_static = ci_p_static, p_dynamic = ci_p_dynamic)
    #   }, error = function(e) {
    #     # Conservative fallback
    #     list(p_static = c(0, 1), p_dynamic = c(0, 1))
    #   }, warning = function(w) {
    #     # Conservative fallback
    #     list(p_static = c(0, 1), p_dynamic = c(0, 1))
    #   })
    # })

    
    # optim_prev<-optim(runif(1),fn=prev_likelihood_real,f_kde=f_kde,n_pools=b_pools,n_sample_size=sample_size_i,
    # observed_ct_vec=observed_ct_vec,neg_lik=TRUE,method = "Brent",lower=0,upper=1)$par

    
      # Extract values 
     optim_prev_esti <- rbind(
      optim_prev_esti,
      data.frame(
        week      = i,
        optim_prev= bbmle::coef(fit)[["ini_p"]],
        optim_p_static   =p_static_hat,
        optim_p_dynamic=p_dynamic_hat,
        # optim_p_static   = bbmle::coef(fit_a)[["p_static"]],
        # optim_p_dynamic   = bbmle::coef(fit_a)[["p_dynamic"]],
        p_q1      = as.numeric(confint1[1]),  # lower bound for ini_p
        p_q2      =  as.numeric(confint1[2]),  # upper bound for ini_p
        static_q1      = ci_p_static[1],
        static_q2      = ci_p_static[2],
        dynamic_q1      = ci_p_dynamic[1],
        dynamic_q2      = ci_p_dynamic[2]
      )
    )
    
    
  }else{
    optim_prev_esti<-rbind(optim_prev_esti,data.frame("week"=i,"optim_prev"=NA,"optim_p_static"=NA,
                                                      "optim_p_dynamic"=NA,
                                                      "p_q1"=NA,"p_q2"=NA,
                                                      "static_q1"=NA,"static_q2"=NA,
                                                       "dynamic_q1"=NA,"dynamic_q2"=NA
                                                      ))
  }
}


optim_prev_esti$total_prev<-optim_prev_esti$optim_p_static+optim_prev_esti$optim_p_dynamic
optim_prev_esti$prop_static<-optim_prev_esti$optim_p_static/optim_prev_esti$total_prev
optim_prev_esti$prop_dynamic<-1-optim_prev_esti$prop_static
optim_prev_esti$total_q1<-optim_prev_esti$static_q1+optim_prev_esti$dynamic_q1
optim_prev_esti$total_q2<-optim_prev_esti$static_q2+optim_prev_esti$dynamic_q2
optim_prev_esti$prop_s_q1<-optim_prev_esti$static_q1/optim_prev_esti$total_q1
optim_prev_esti$prop_s_q2<-optim_prev_esti$static_q2/optim_prev_esti$total_q2

ggplot(optim_prev_esti)+
  geom_line(aes(x=week,y=optim_prev),color="black")+
  geom_line(aes(x=week,y=total_prev),color="red")

p1<-ggplot(optim_prev_esti)+
  geom_ribbon(aes(x=week,ymin=static_q1,ymax = static_q2),fill="#d2b48c",alpha=0.7)+
  geom_point(aes(x=week,y=optim_p_static),color="black")+
ylim(0,1)


p2<-ggplot(optim_prev_esti)+
  geom_ribbon(aes(x=week,ymin=dynamic_q1,ymax = dynamic_q2),fill="#d2b48c",alpha=0.7)+
  geom_point(aes(x=week,y=optim_p_dynamic),color="black")+
  ylim(0,1)

p1|p2

ggplot(optim_prev_esti)+
 geom_point(aes(x=week,y=prop_dynamic),color="darkblue")+
  geom_point(aes(x=week,y=prop_static),color="green")+
geom_hline(yintercept = .5, color = "red", linetype = "dashed", size = 1)+
  ylab("Estimated proprtion of static and dynamic")
  #geom_point(aes(x=week,y=prop_dynamic),color="red")+
  # geom_ribbon(data=optim_prev_esti,aes(x=week,ymin=prop_s_q1,ymax=prop_s_q2),
  #             fill="#d2b48c",alpha=0.7)+
 # ylim(0,2)
  


prev_all_dat<-subset(all_data, all_data$year==yr_to_run & all_data$sim_number==sim_to_run)

prev_data <- prev_all_dat %>% mutate(positive=ct_value<ct_threshold) %>% group_by(current_time) %>%
  summarize(prev=sum(positive)/n(),N=n(),prev_true=sum(inf_state)/n())


all_prev<-data.frame("week"=1:52,"actual_prev"=prev_data$prev_true,"estimated_prev"=optim_prev_esti$optim_prev)
m_all_prev<-melt(all_prev,id="week")
col_names<-c("Actual prevalence", "Estimated prevalence \n (with pooled Ct values)")
levels(m_all_prev$variable)<-col_names


p_est_prev<-ggplot()+
  geom_ribbon(data=optim_prev_esti,aes(x=week,ymin=p_q1,ymax=p_q2),fill="#d2b48c",alpha=0.7)+
  geom_line(data=m_all_prev,aes(x=week,y=value,color=variable),linewidth=1,alpha=.9)+
  scale_color_manual(values = c("#009ff7","#9d6723","#005841"),labels=col_names)+
  # scale_color_manual(values=c("black","#008e7c"))+
  #geom_point(data=act_prev,aes(x=week,y=actual_prev),color="#7c2d1d")+
  # geom_point(data=opt_prev,aes(x=week,y=estimated_prev),color="#9d6723",size=2)+
  default_theme+
  # scale_color_discrete_qualitative(palette = "Dark3")+
  ylim(0,0.1)+
  xlab("Time (week)")+
  ylab("Prevalence")+
  theme(legend.position = "bottom")+
  ggtitle("(C) Estimated prevalence \n using pooled Ct values")
p_est_prev

p_proportion<-ggplot()+
  geom_ribbon(data=optim_prev_esti,aes(x=week,ymin=prop_s_q1,ymax=prop_s_q2),fill="#d2b48c",alpha=0.7)+
  geom_point(data=optim_prev_esti,aes(x=week,y=prop_static),linewidth=1,alpha=.9)+
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") +
 # geom_hline(yintercept = 0.4, linetype = "dotted", color = "black") +
#  geom_hline(yintercept = 0.6, linetype = "dotted", color = "black") +
  #scale_color_manual(values = c("#009ff7","#9d6723","#005841"),labels=col_names)+
  # scale_color_manual(values=c("black","#008e7c"))+
  #geom_point(data=act_prev,aes(x=week,y=actual_prev),color="#7c2d1d")+
  # geom_point(data=opt_prev,aes(x=week,y=estimated_prev),color="#9d6723",size=2)+
  default_theme+
  # scale_color_discrete_qualitative(palette = "Dark3")+
 ylim(0,1)+
  xlab("Time (week)")+
  ylab("Proportion")+
  theme(legend.position = "bottom")
#  ggtitle("(C) Estimated prevalence \n using pooled Ct values")

p_proportion

#cdc package 
pv_cdc<-data.frame(matrix(NA,nrow = 0,ncol = 4))
colnames(pv_cdc)<-c("week","cdc_prev","q1","q2")
inf_samples<-NULL

#calculate the prevalence from the cdc package
for(i in 1:52){
  cts_at_i<-subset(cts_per_poolss,cts_per_poolss$week==i)
  inf_states<-(cts_at_i$pooled_ct_per_pool<40)*1
  pol_sz<-cts_at_i$sample_size
  inf_samples<-rbind(inf_samples,data.frame("week"=i,"postives"=sum(inf_states),"negatives"=length(inf_states)-sum(inf_states)))
  if(sum(pol_sz)>0){
    pv<-pooledBin(inf_states,pol_sz,ci.method ="lrt" )
    pv_cdc<-rbind(pv_cdc,data.frame("week"=i,"cdc_prev"=pv$P,"q1"=pv$Lower,"q2"=pv$Upper))
  }
  else{
    pv_cdc<-rbind(pv_cdc,data.frame("week"=i,"cdc_prev"=NA,"q1"=NA,"q2"=NA))
  }
}


all_prev<-data.frame("week"=1:52,"actual_prev"=prev_data$prev_true,"estimated_prev"=pv_cdc$cdc_prev)
m_all_prev<-melt(all_prev,id="week")
col_names<-c("Actual prevalence", "Estimated prevalence \n (with binarised data)")
levels(m_all_prev$variable)<-col_names

p_est_prev_cdc<-ggplot()+
  geom_ribbon(data=pv_cdc,aes(x=week,ymin=q1,ymax=q2),fill="#005841",alpha=0.4)+
  geom_line(data=m_all_prev,aes(x=week,y=value,color=variable),linewidth=1,alpha=.9)+
  scale_color_manual(values = c("#009ff7","#005841"),labels=col_names)+
  # scale_color_manual(values=c("black","#008e7c"))+
  #geom_point(data=act_prev,aes(x=week,y=actual_prev),color="#7c2d1d")+
  # geom_point(data=opt_prev,aes(x=week,y=estimated_prev),color="#9d6723",size=2)+
  default_theme+
  # scale_color_discrete_qualitative(palette = "Dark3")+
  ylim(0,0.1)+
  xlab("Time (week)")+
  ylab("Prevalence")+
  theme(legend.position = "bottom")+
  ggtitle("(D) Estimated prevalence \n with binarised data")
p_est_prev_cdc


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
  ggtitle("(A) Pooled Ct values")

ct_per_pools

m_inf_samples<-melt(inf_samples,id="week")
m_inf_samples=subset(m_inf_samples,week %in% 18:42)

pos_pl<-ggplot(data=m_inf_samples,aes(as.factor(week),y=value,fill=variable))+
  geom_bar(stat="identity")+
  default_theme+
  xlab("Time (weeks)")+
  ylab("Count")+
  theme(legend.position = "bottom")+
  scale_fill_manual(values  =c("#903920","#759228"),,labels= c("Positive", "Negative"))+
  ggtitle("(B) Binarised data")
pos_pl



(ct_per_pools +pos_pl)/ (p_est_prev+p_est_prev_cdc)

#ggsave("figures/figure_5.png",last_plot(),height = 8,width = 12)
