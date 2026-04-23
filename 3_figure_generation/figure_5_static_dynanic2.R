#load the default theme for figures:
theme<-source("functions/default_theme.R")
default_theme<-default_theme()

#if PooledInfRate is not installed: 

#library(devtools)
#devtools::install("/Users/palahakoon/Documents/GitHub/PooledInfRate", build_vignettes = FALSE)

source("functions/prev_likelihood_real_static_dynamic5_vec_reparametrise2.R")
source("functions/prev_likelihood_real.R")
#load all the data:
all_data <- readRDS("synthetic_data/all_data.rds")
#load the data: pooled ct values 
cts_per_pool= readRDS("pre_calculations/pooled_cts/cts_per_pool_demo.rds") #they are only data from 2:4 years for the 4 sims, 

ct_threshold=40


#also doubke check if the sample size is greater than 50:
cts_per_pool$sample_size[cts_per_pool$sample_size>50]=50



#choose the sim number and the year to run:
yr_to_run=4
sim_to_run=2
cts_per_poolss<-subset(cts_per_pool,(cts_per_pool$sim_number==sim_to_run) & (cts_per_pool$year==yr_to_run))


#read all the fkdes 
all_f_kde=NULL
all_f_kde=sapply(1:50,function(x) readRDS(file = paste0("pre_calculations/kde_calcs/f_kde_",x,".rds")))

#joint KDEs taht were precalculated. See KDE calculations to see how they were pre-caluculated
all_f_kde_joint=NULL
all_f_kde_joint=lapply(1:50,function(x) readRDS(file = paste0("pre_calculations/joint_kde_calc/joint_fkde_",x,".rds")))


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
        if(is.na(ci[2])) ci[2] <- 1
        ci
      }, error = function(e) {
        c(0, 1)
      }, warning = function(w) {
        c(0, 1)
      })
    )
    
    #to calculate the the dynamc vs static proposrions, do it independently: 
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
    
    # Extract values 
    optim_prev_esti <- rbind(
      optim_prev_esti,
      data.frame(
        week      = i,
        optim_prev= bbmle::coef(fit)[["ini_p"]],
        prev_q1=confint1[1],
        prev_q2=confint1[2],
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
    optim_prev_esti<-rbind(optim_prev_esti,data.frame("week"=i,"optim_prev"=NA,
                                                      "prev_q1"=NA,
                                                      "prev_q2"=NA,
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
    ))
  }
}


#to make the. new Figure 5, run these codes: 
optim_prev_esti$total_prev<-optim_prev_esti$optim_p_static+optim_prev_esti$optim_p_dynamic
optim_prev_esti$prop_static<-optim_prev_esti$optim_p_static/optim_prev_esti$total_prev
optim_prev_esti$prop_dynamic<-1-optim_prev_esti$prop_static



prev_all_dat<-subset(all_data, all_data$year==yr_to_run & all_data$sim_number==sim_to_run)



prev_data <- prev_all_dat %>%
  mutate(
    positive = ct_value < ct_threshold,
    dynamic_inf = (inf_state == 1 & ct_type_method == 3)
  ) %>%
  group_by(current_time) %>%
  summarise(
    prev = sum(positive) / n(),
    N = n(),
    prev_true = sum(inf_state) / n(),                      # original true prevalence
    prev_true_dynamic = sum(dynamic_inf) / n()              # NEW dynamic true prevalence
  )


all_prev<-data.frame("week"=1:52,"actual_prev"=prev_data$prev_true,"prev_true_dynamic"=prev_data$prev_true_dynamic)
all_prev<-cbind(all_prev,optim_prev_esti)
all_prev<-all_prev[,-4]



c_dat<-data.frame("week"=all_prev$week,"actual_prev"=all_prev$actual_prev,"optim_prev"=all_prev$optim_prev)
m_c_dat<-melt(c_dat,id="week")

col_names<-c("Actual prevalence", "Estimated prevalence \n (with pooled Ct data)")

# ensure variable is a factor in the desired order and drop unused levels
m_c_dat$variable <- factor(m_c_dat$variable, levels = unique(m_c_dat$variable))
m_c_dat$variable <- droplevels(m_c_dat$variable)

p_D <- ggplot() +
  geom_ribbon(
    data = all_prev,
    aes(x = week, ymin = prev_q1, ymax = prev_q2),
    fill = "#d2b48c",
    alpha = 0.7
  ) +
  geom_line(
    data = m_c_dat,
    aes(
      x = week,
      y = value,
      color = variable,
      linetype = variable         # <- map linetype here
    ),
    linewidth = 1,
    alpha = .9
  ) +
  scale_color_manual(
    name   = NULL,               # match linetype name
    values = c("#009ff7", "#9d6723"),
    labels = col_names
  ) +
  scale_linetype_manual(
    name   = NULL,             
    values = c("dotdash", "solid"),
    labels = col_names
  ) +
  default_theme +
  ylim(0, 0.1) +
  xlab("Time (week)") +
  ylab("Prevalence") +
  theme(legend.position = "bottom") +
  ggtitle("(D) Estimated prevalence \n using pooled Ct values")

p_D


e_dat<-data.frame("week"=all_prev$week,"actual_prev"=all_prev$actual_prev,"prev_true_dynamic"=all_prev$prev_true_dynamic,
                  "optim_p_dynamic"=all_prev$optim_p_dynamic)
m_e_dat<-melt(e_dat,id="week")
col_names<-c("Actual prevalence", "Actual productive\n infection prevalence","Estimated productive\n infection prevalence \n (with pooled Ct data)")


m_e_dat$variable <- factor(m_e_dat$variable, levels = unique(m_e_dat$variable))

p_F <- ggplot() +
  geom_ribbon(
    data = all_prev,
    aes(x = week, ymin = dynamic_q1, ymax = dynamic_q2),
    fill = "#d8bfd8", alpha = 0.7
  ) +
  geom_line(
    data = m_e_dat,
    aes(x = week, y = value, color = variable, linetype = variable),
    linewidth = 1, alpha = .9
  ) +
  scale_color_manual(
    name   = NULL,                       
    values = c("#009ff7", "#660049", "#b23d8c"),
    labels = col_names
  ) +
  scale_linetype_manual(
    name   = NULL,                   
    values = c("dotdash", "solid", "solid"),
    labels = col_names
  ) +
  default_theme +
  ylim(0, 0.1) +
  xlab("Time (week)") +
  ylab("Prevalence") +
  theme(legend.position = "bottom") +
  ggtitle("(F) Estimated productive infection prevalence \n using pooled Ct values")

p_F



# 1) Long data for the line
new_dat <- tibble::tibble(
  week         = all_prev$week,
  prop_static  = all_prev$prop_static,
  prop_dynamic = all_prev$prop_dynamic
)

line_long <- new_dat %>%
  pivot_longer(cols = c(prop_static, prop_dynamic),
               names_to = "method",
               values_to = "prop") %>%
  mutate(method = ifelse(method == "prop_static", "Static", "Dynamic"))

# 2) Quantiles: keep one row per week-method with q1/q2
qs_long <- tibble::tibble(
  week            = all_prev$week,
  prop_static_q1  = all_prev$prop_static_q1,
  prop_static_q2  = all_prev$prop_static_q2,
  prop_dynamic_q1 = all_prev$prop_dynamic_q1,
  prop_dynamic_q2 = all_prev$prop_dynamic_q2
) %>%
  pivot_longer(cols = c(prop_static_q1, prop_static_q2, prop_dynamic_q1, prop_dynamic_q2),
               names_to = c("method", ".value"),
               names_pattern = "prop_(static|dynamic)_(q[12])") %>%
  # After names_pattern, columns are q1 and q2; method is 'static'/'dynamic'
  mutate(method = ifelse(method == "static", "Static", "Dynamic"))

# 3) Join lines with their matching ribbons
plot_dat <- line_long %>%
  inner_join(qs_long, by = c("week", "method"))
# plot_dat has: week, method, prop (line), q1 (ymin), q2 (ymax)


plot_dat$method <- factor(plot_dat$method,
                          levels = c("Static", "Dynamic"),
                          labels = c("Non productive", "Productive"))



p_C <- ggplot(subset(plot_dat,plot_dat$method=="Productive"), aes(x = week)) +
  geom_ribbon(aes(ymin = q1, ymax = q2, fill = method), alpha = 0.25, colour = NA) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black", linewidth = 0.7) +
  geom_line(aes(y = prop, colour = method), linewidth = 1) +
 # facet_wrap(~ method, ncol = 2, scales = "free_y") +
  labs(x = "Time (week)", y = "Proportion") +
  default_theme +
  theme(legend.position = "none") +
  ggtitle("(C) Proportion of productive \n infections in the population")

p_C




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


all_prev<-data.frame("week"=1:52,"actual_prev"=all_prev$actual_prev,"estimated_prev"=pv_cdc$cdc_prev)
m_all_prev<-melt(all_prev,id="week")
m_all_prev <- all_prev %>%
  pivot_longer(cols = -week, names_to = "variable", values_to = "value")


# 2) pivot longer
m_all_prev <- all_prev %>%
  pivot_longer(cols = -week, names_to = "variable", values_to = "value")

# 3) define nice labels and force factor ordering
col_names <- c("Actual prevalence", "Estimated prevalence \n (with binarised data)")
m_all_prev$variable <- factor(m_all_prev$variable,
                              levels = c("actual_prev", "estimated_prev"),
                              labels = col_names)

p_E <- ggplot() +
  geom_ribbon(data = pv_cdc, aes(x = week, ymin = q1, ymax = q2),
              fill = "#005841", alpha = 0.4) +
  geom_line(data = m_all_prev,
            aes(x = week, y = value, color = variable, linetype = variable),
            linewidth = 1, alpha = .9) +
  scale_color_manual(name = NULL,
                     values = c("#009ff7", "#005841"),
                     labels = col_names) +
  scale_linetype_manual(name = NULL,
                        values = c("dotdash", "solid"),
                        labels = col_names) +
  default_theme +
  ylim(0, 0.1) +
  xlab("Time (week)") +
  ylab("Prevalence") +
  theme(legend.position = "bottom") +
  ggtitle("(E) Estimated prevalence \n with binarised data")

p_E

p_A<-ggplot(data=subset(cts_per_pool,cts_per_pool$year==yr_to_run & cts_per_pool$sim_number==sim_to_run),aes(x=(week),y=pooled_ct_per_pool))+
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

p_A


m_inf_samples<-melt(inf_samples,id="week")
m_inf_samples=subset(m_inf_samples,week %in% 18:42)

p_B<-ggplot(data=m_inf_samples,aes(as.factor(week),y=value,fill=variable))+
  geom_bar(stat="identity")+
  default_theme+
  xlab("Time (weeks)")+
  ylab("Count")+
  theme(legend.position = "bottom")+
  scale_fill_manual(values  =c("#903920","#759228"),,labels= c("Positive", "Negative"))+
  ggtitle("(B) Binarised data")
p_B

(p_A+p_B+p_C)/ (p_D+p_E+p_F)


ggsave("figures/figure_5_static_dynamic.png",last_plot(),height = 10,width = 20)




###cpmparison of the two methods; 

