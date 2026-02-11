
#load the default theme for figures:
  theme<-source("functions/default_theme.R")
  default_theme<-default_theme()
  
#write a script to check when the pool sizes are long  prevelance of the two are high:

#generate a known prevlanace:
 #take pool size to be 50 
pool_sizes=50 
n_pools<-c(10,30,50,100)
      
#get the proprtion sample 
#to have re
p_static=c(0.001,0.01,.1)
p_dynamic=c(0.001,0.01,.1)
#p_dynamic=.2 


# p_combi<-expand.grid(p_static,p_dynamic,pool_sizes,n_pools)
# 
# 

grid <- expand.grid(p_static, p_dynamic, pool_sizes, n_pools)
p_combi <- crossing(rep = 1:100, grid)
n_prev<-nrow(p_combi)
n_prev

ct_threshold=40

#generate random observed Ct vectors 
all_data<-readRDS("synthetic_data/all_data.rds")


#get the ct values that are static and and dynamic 


static_moz<-subset(all_data,ct_type_method==2)
dynamic_moz<-subset(all_data,ct_type_method==3)

#KDE estimation:
positive_ctss_static<-static_moz %>% filter(ct_value < ct_threshold) %>% pull(ct_value)
positive_ctss_dynamic<-dynamic_moz %>% filter(ct_value < ct_threshold) %>% pull(ct_value)

#transform to viral loads:
ct_to_vl<-function(ct,intercept,slope){ #ct=m*log10(vl)+c0
  vl<-10^((ct-intercept)/slope)
  vl   
}
intercept=36.9
slope= -2.7 


viral_loads_vec_static=ct_to_vl(positive_ctss_static,intercept,slope)
viral_loads_vec_dynamic=ct_to_vl(positive_ctss_dynamic,intercept,slope)


ct_vec_obs<-function(prev,pool_sizes,n_pools,viral_loads_vec_static,viral_loads_vec_dynamic){
  
  vl_to_ct<-function(slope,intercept,vl){ #vl in 10^
    ct1<-slope*log10(vl) +intercept
    ct<-min(ct1,40)
    ct
  }
  
  intercept=36.9
  slope= -2.7 
  
  all_pool_cts<-NULL
  
  for (i in 1:n_pools) {
    n_positives_static<-rbinom(1,pool_sizes,as.numeric(prev[1]))
    n_positives_dynamic<-rbinom(1,pool_sizes,as.numeric(prev[2]))
    pos_vls_stat<- sample(viral_loads_vec_static,n_positives_static)
    pos_vls_dyna <- sample(viral_loads_vec_dynamic,n_positives_dynamic)
    pos_vls<-c(pos_vls_stat,pos_vls_dyna)
    neg_vls <- rep(0, pool_sizes-n_positives_static-n_positives_dynamic)
    pool_vl<- mean(c(pos_vls,neg_vls))
    pool_ct<-vl_to_ct(slope,intercept,pool_vl)
    all_pool_cts[i]<-pool_ct
  }

  return(all_pool_cts)
}


all_pools<-NULL #each list has the observed pooled Ct values given the number of poools 
#for each combination of prevalance,

for (i in 1:n_prev) {
    prv<-as.numeric(p_combi[i,2:3])
    pool_s<-as.numeric(p_combi[i,4])
    n_pool<-as.numeric(p_combi[i,5])
    pool_ct<-ct_vec_obs(prv,pool_sizes = pool_s,n_pools = n_pool,viral_loads_vec_static,viral_loads_vec_dynamic)
    all_pools[[i]]<-pool_ct
  
}

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

#rows<-prevelance, columns, each pooled Ct for each pool 
#<-as.data.frame(all_pools)

#calculate the prevalances: 
all_f_kde_joint=NULL
all_f_kde_joint=lapply(1:50,function(x) readRDS(file = paste0("pre_calculations/joint_kde_calc/joint_fkde_",x,".rds")))

# 
# all_f_kde_static=NULL
# all_f_kde_static=lapply(1:50,function(x) readRDS(file = paste0("pre_calculations/kde_calcs_dynamic_static/f_kde_static_",x,".rds")))
# 
# all_f_kde_dynamic=NULL
# all_f_kde_dynamic=lapply(1:50,function(x) readRDS(file = paste0("pre_calculations/kde_calcs_dynamic_static/f_kde_dynamic_",x,".rds")))

#run the prevalence estimation: 
optim_prev_esti=NULL
S <- 5000



is_pd_chol <- function(Sigma) {
  # expects a symmetric matrix
  out <- try(chol(Sigma), silent = TRUE)
  !inherits(out, "try-error")
}

make_pd <- function(Sigma) {
  # 1) symmetrize
  S <- (Sigma + t(Sigma)) / 2
  
  # 2) return as-is if already PD
  if (is_pd_chol(S)) return(S)
  
  # 3) otherwise, project to nearest PD
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for nearPD(). Please install it.")
  }
  S_pd <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
  return(S_pd)
}




nn<-n_prev #n_prev
for (i in 1612:nn) {
  b_pools<-as.numeric(p_combi[i,5])
  sample_size_i<-as.numeric(rep(p_combi[i,4],b_pools))
  #observed_ct_vec<-as.numeric(all_pools[i,])
  observed_ct_vec<-all_pools[[i]]
    #to calculate the the dynamc vs static proposrions, do it independently: 
    # f_kde_static=NULL
    # f_kde_static=all_f_kde_static[sample_size_i] #get the kde's sample sizes 
    # f_kde_dynamic=NULL
    # f_kde_dynamic=all_f_kde_dynamic[sample_size_i] #get the kde's sample sizes 
    
    # 
    # fit_a <- bbmle::mle2(prev_likelihood_real_static_dynamic5_vec,start=list(p_static=runif(1),
    #                                                                          p_dynamic=runif(1)),
    #                      data=list(f_kde_static=f_kde_static,
    #                                f_kde_dynamic=f_kde_dynamic,
    #                                n_pools=b_pools,
    #                                pool_sizes=sample_size_i,
    #                                observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
    #                      method="L-BFGS-B",
    #                      lower=c(0,0),upper=c(1,1))
    # 
    
    fit_a <- bbmle::mle2(prev_likelihood_real_static_dynamic5_vec_reparametrise2,
                         start=list(a = qlogis(0.005), b = qlogis(0.005)),
                         data=list(f_kde_joint=all_f_kde_joint,all_combos=all_combos,
                                   n_pools=b_pools,
                                   pool_sizes=sample_size_i,
                                   observed_ct_vec=observed_ct_vec,neg_lik=TRUE),
                         method="L-BFGS-B",
                         lower=c(-50,-50),upper=c(50,50))
    a_hat <- coef(fit_a)[["a"]]
    b_hat <- coef(fit_a)[["b"]]
    p_static_hat <- plogis(a_hat)
    p_dynamic_hat <- plogis(b_hat) * (1 - p_static_hat)
    
    cov_mat<-make_pd(vcov(fit_a)) 
    boots <- MASS::mvrnorm(S, mu = coef(fit_a), Sigma = cov_mat)
    a_sim <- boots[,1]
    b_sim <- boots[,2]
    
    p_static_sim  <- plogis(a_sim)
    p_dynamic_sim <- plogis(b_sim) * (1 - p_static_sim)
    #also calculate the proportions:
    p_dynamic_sim_prop<-p_dynamic_sim/(p_dynamic_sim+p_static_sim)
    p_static_sim_prop<-p_static_sim/(p_dynamic_sim+p_static_sim)
    
    #confidence intervals 
    ci_p_static  <- quantile(p_static_sim,  probs = c(.025, .975))
    ci_p_dynamic <- quantile(p_dynamic_sim, probs = c(.025, .975))
    
    ci_prop_static<-quantile(p_static_sim_prop,  probs = c(.025, .975),na.rm=T)
    ci_prop_dynamic<-quantile(p_dynamic_sim_prop,  probs = c(.025, .975),na.rm=T)
    # suppressMessages({
    #   confint_list <- tryCatch({
    #     # Profile likelihood CIs for both parameters
    #     ci_all <- bbmle::confint(fit_a, method = "quad")
    #     #or use this for a a more accurate but reduced number of steps
    #     #using the profilinh method: 
    #     # ci_all <- bbmle::confint(fit_a, method = "profile", nsteps = 50)
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
        optim_p_static   = p_static_hat,
        optim_p_dynamic   = p_dynamic_hat,
        static_q1      = ci_p_static[1],
        static_q2      = ci_p_static[2],
        dynamic_q1      = ci_p_dynamic[1],
        dynamic_q2      = ci_p_dynamic[2],
        prop_dynamic_q1=ci_prop_dynamic[1],
        prop_dynamic_q2=ci_prop_dynamic[2],
        prop_static_q1=ci_prop_static[1],
        prop_static_q2=ci_prop_static[2]
        # optim_p_static   = bbmle::coef(fit_a)[["p_static"]],
        # optim_p_dynamic   = bbmle::coef(fit_a)[["p_dynamic"]]
        # static_q1      = confint_list$p_static[1],
        # static_q2      = confint_list$p_static[2],  
        # dynamic_q1      = confint_list$p_dynamic[1],
        # dynamic_q2      = confint_list$p_dynamic[2]
      )
    )
    
    
 
}


saveRDS(optim_prev_esti, "optim_prev_esti_5.rds")

# 
# 
# 
# # pre-allocate if not already
# if (!exists("optim_prev_esti")) {
#   optim_prev_esti <- data.frame()
# }
# 
# for (i in 1:nn) {
#   b_pools <- p_combi[i, 4]
#   sample_size_i <- rep(p_combi[i, 3], b_pools)
#   observed_ct_vec <- all_pools[[i]]
#   
#   fit_a <- bbmle::mle2(
#     prev_likelihood_real_static_dynamic5_vec_reparametrise2,
#     start = list(a = qlogis(0.005), b = qlogis(0.005)),
#     data = list(
#       f_kde_joint = all_f_kde_joint,
#       all_combos   = all_combos,
#       n_pools      = b_pools,
#       pool_sizes   = sample_size_i,
#       observed_ct_vec = observed_ct_vec,
#       neg_lik      = TRUE
#     ),
#     method = "L-BFGS-B",
#     lower = c(-50, -50),
#     upper = c(50, 50)
#   )
#   
#   # Point estimates on parameter and probability scales
#   coef_est <- coef(fit_a)
#   a_hat <- coef_est[["a"]]
#   b_hat <- coef_est[["b"]]
#   
#   p_static_hat  <- plogis(a_hat)
#   p_dynamic_hat <- plogis(b_hat) * (1 - p_static_hat)
#   
#   # Obtain quadratic approximation CIs for parameters a and b
#   conf_res <- tryCatch({
#     bbmle::confint(fit_a, method = "quad")
#   }, error = function(e) {
#     warning(sprintf("confint(method='quad') failed at i=%s: %s", i, conditionMessage(e)))
#     NULL
#   }, warning = function(w) {
#     # still return the object if produced but warn
#     tryCatch(bbmle::confint(fit_a, method = "quad"), error = function(e) NULL)
#   })
#   
#   if (is.null(conf_res) || !all(c("a", "b") %in% rownames(conf_res))) {
#     # conservative fallback
#     a_ci <- c(-Inf, Inf)
#     b_ci <- c(-Inf, Inf)
#     
#     ci_p_static  <- c(0, 1)
#     ci_p_dynamic <- c(0, 1)
#     ci_prop_static  <- c(0, 1)
#     ci_prop_dynamic <- c(0, 1)
#   } else {
#     # conf_res rows named "a" and "b": extract their (lower, upper)
#     a_ci <- as.numeric(conf_res["a", ])
#     b_ci <- as.numeric(conf_res["b", ])
#     
#     # transform parameter CIs to probability scale
#     # NOTE: this applies the transformation to each parameter independently
#     # p_static CI:
#     ci_p_static <- plogis(a_ci) # length-2 numeric vector: lower, upper
#     
#     # p_dynamic depends on both a and b. A direct but approximate approach:
#     # transform b-ci and combine with p_static endpoints. We produce two candidate
#     # endpoints by pairing b lower/upper with a lower/upper in a conservative manner.
#     # Here we compute the dynamic probability at the four corners and take min/max:
#     a_vals <- a_ci
#     b_vals <- b_ci
#     # corners: (a_lo,b_lo), (a_lo,b_hi), (a_hi,b_lo), (a_hi,b_hi)
#     corner_p_dynamic <- vapply(
#       list(
#         c(a_vals[1], b_vals[1]),
#         c(a_vals[1], b_vals[2]),
#         c(a_vals[2], b_vals[1]),
#         c(a_vals[2], b_vals[2])
#       ),
#       FUN.VALUE = numeric(1),
#       FUN = function(ab) {
#         a_val <- ab[1]; b_val <- ab[2]
#         p_s <- plogis(a_val)
#         p_d <- plogis(b_val) * (1 - p_s)
#         return(p_d)
#       }
#     )
#     ci_p_dynamic <- c(min(corner_p_dynamic, na.rm = TRUE), max(corner_p_dynamic, na.rm = TRUE))
#     
#     # Derived proportions: dynamic / (dynamic + static)
#     # compute the proportions at the same four corners and take min/max
#     corner_prop_dynamic <- vapply(
#       list(
#         c(a_vals[1], b_vals[1]),
#         c(a_vals[1], b_vals[2]),
#         c(a_vals[2], b_vals[1]),
#         c(a_vals[2], b_vals[2])
#       ),
#       FUN.VALUE = numeric(1),
#       FUN = function(ab) {
#         a_val <- ab[1]; b_val <- ab[2]
#         p_s <- plogis(a_val)
#         p_d <- plogis(b_val) * (1 - p_s)
#         denom <- p_d + p_s
#         if (denom <= 0 || !is.finite(denom)) return(NA_real_)
#         return(p_d / denom)
#       }
#     )
#     corner_prop_static <- vapply(
#       list(
#         c(a_vals[1], b_vals[1]),
#         c(a_vals[1], b_vals[2]),
#         c(a_vals[2], b_vals[1]),
#         c(a_vals[2], b_vals[2])
#       ),
#       FUN.VALUE = numeric(1),
#       FUN = function(ab) {
#         a_val <- ab[1]; b_val <- ab[2]
#         p_s <- plogis(a_val)
#         p_d <- plogis(b_val) * (1 - p_s)
#         denom <- p_d + p_s
#         if (denom <= 0 || !is.finite(denom)) return(NA_real_)
#         return(p_s / denom)
#       }
#     )
#     
#     ci_prop_dynamic <- c(min(corner_prop_dynamic, na.rm = TRUE), max(corner_prop_dynamic, na.rm = TRUE))
#     ci_prop_static  <- c(min(corner_prop_static,  na.rm = TRUE), max(corner_prop_static,  na.rm = TRUE))
#     
#     # In case any NA forced min/max to be Inf/-Inf:
#     if (!is.finite(ci_p_dynamic[1])) ci_p_dynamic[1] <- 0
#     if (!is.finite(ci_p_dynamic[2])) ci_p_dynamic[2] <- 1
#     if (!is.finite(ci_prop_dynamic[1])) ci_prop_dynamic[1] <- 0
#     if (!is.finite(ci_prop_dynamic[2])) ci_prop_dynamic[2] <- 1
#     if (!is.finite(ci_prop_static[1])) ci_prop_static[1] <- 0
#     if (!is.finite(ci_prop_static[2])) ci_prop_static[2] <- 1
#   }
#   
#   # Append to results (fixed mapping for prop_static_q2)
#   optim_prev_esti <- rbind(
#     optim_prev_esti,
#     data.frame(
#       week            = i,
#       optim_p_static  = p_static_hat,
#       optim_p_dynamic = p_dynamic_hat,
#       static_q1       = ci_p_static[1],
#       static_q2       = ci_p_static[2],
#       dynamic_q1      = ci_p_dynamic[1],
#       dynamic_q2      = ci_p_dynamic[2],
#       prop_dynamic_q1 = ci_prop_dynamic[1],
#       prop_dynamic_q2 = ci_prop_dynamic[2],
#       prop_static_q1  = ci_prop_static[1],
#       prop_static_q2  = ci_prop_static[2]
#     )
#   )
# }
# 
# 
# 


#saveRDS(optim_prev_esti, "optim_prev_esti_4.rds")



#use optim_prev_esti_2.rds to generate figures 1 and 2
optim_prev_esti<-readRDS("optim_prev_estims.rds")


optim_prev_esti$true_static_p<-as.vector(p_combi[,2])
optim_prev_esti$true_dynamic_p<-as.vector(p_combi[,3])
optim_prev_esti$true_dynamic_prop<-optim_prev_esti$true_dynamic_p/(optim_prev_esti$true_static_p+optim_prev_esti$true_dynamic_p)
optim_prev_esti$estimated_dynamic_prop<-optim_prev_esti$optim_p_dynamic/(optim_prev_esti$optim_p_dynamic+optim_prev_esti$optim_p_static)
optim_prev_esti$propr_error_dynamic<-optim_prev_esti$true_dynamic_prop-optim_prev_esti$estimated_dynamic_prop

optim_prev_esti$n_pools<-p_combi[,5]
optim_prev_esti$pool_size<-p_combi[,4]

#remove 30, 50 samples as they are repeated from the above:

# optim_prev_esti <- subset(optim_prev_esti, !n_pools %in% c(30, 50))



# optim_prev_esti<-rbind(optim_prev_esti,optim_prev_esti_2)

# saveRDS(optim_prev_esti, "optim_prev_esti_final.rds")

lvs=paste0("Pool size = ",pool_sizes)
optim_prev_esti$pool_size=as.factor(optim_prev_esti$pool_size)
levels(optim_prev_esti$pool_size)=lvs

lv2=paste0("Prevalence = ",p_dynamic)
optim_prev_esti$true_Prevalence=as.factor(optim_prev_esti$true_dynamic_p)
levels(optim_prev_esti$true_Prevalence)=lv2



#plot the estimates with actual data 

new_dat<-subset(optim_prev_esti, pool_size == "Pool size = 50")

levs<-round(unique( new_dat$true_dynamic_prop),2)
new_dat$true_dynamic_prop<-as.factor(new_dat$true_dynamic_prop)
levels(new_dat$true_dynamic_prop)<-levs


p1 <- ggplot(data = new_dat) +
  geom_point(
    aes(
      x = as.factor(true_dynamic_prop),
      y = propr_error_dynamic,
      colour = as.factor(n_pools)
    )
  ) +
  geom_abline(
    slope = 0,
    intercept = 0,
    colour = "black",
    linetype = "dashed"
  ) +
  facet_grid(~ as.factor(true_Prevalence)) +
  labs(
    x = "True dynamic proportion",
    y = "Dynamic prevalence estimation error",
    colour = "Number of pools"     
  ) +
  theme(legend.position = "bottom")

p1

p2 <- ggplot(data = new_dat) +
  geom_errorbar(
    aes(x= as.factor(true_dynamic_prop),
              ymin = dynamic_q1,
              ymax = dynamic_q2,
              colour = as.factor(n_pools)),
    position = position_dodge(width =0.4)
  )+
  geom_point(
    aes(
      x = as.factor(true_dynamic_prop),
      y = optim_p_dynamic,
      colour = as.factor(n_pools)
    ),     position = position_dodge(width = .4)
  ) +
  facet_wrap(~ as.factor(true_Prevalence),scales="free_y") +
  labs(
    x = "True dynamic proportion",
    y = "Estimated prevalence (productive infection) ",
    colour = "Number of pools"     
  ) +
  theme(legend.position = "bottom")

p2




###use optim_prev_esti_5.rds to geretae figure 3 
optim_prev_esti_2<-as.data.frame(readRDS("optim_prev_estims.rds"))

optim_prev_esti_2$true_static_p<-as.vector(p_combi[,2][[1]])
optim_prev_esti_2$true_dynamic_p<-as.vector(p_combi[,3][[1]])
optim_prev_esti_2$true_dynamic_prop<-optim_prev_esti_2$true_dynamic_p/(optim_prev_esti_2$true_static_p+optim_prev_esti_2$true_dynamic_p)
optim_prev_esti_2$estimated_dynamic_prop<-optim_prev_esti_2$optim_p_dynamic/(optim_prev_esti_2$optim_p_dynamic+optim_prev_esti_2$optim_p_static)
optim_prev_esti_2$propr_error_dynamic<-optim_prev_esti_2$true_dynamic_prop-optim_prev_esti_2$estimated_dynamic_prop

optim_prev_esti_2$n_pools<-as.vector(p_combi[,5][[1]])
optim_prev_esti_2$pool_size<-as.vector(p_combi[,4][[1]])


lvs=paste0("Pool size = ",pool_sizes)
optim_prev_esti_2$pool_size=as.factor(optim_prev_esti_2$pool_size)
levels(optim_prev_esti_2$pool_size)=lvs
 
lv2=paste0("Prevalence = ",p_dynamic)
optim_prev_esti_2$true_Prevalence=as.factor(optim_prev_esti_2$true_dynamic_p)
levels(optim_prev_esti_2$true_Prevalence)=lv2



#plot the estimates with actual data 

new_dat<-subset(optim_prev_esti_2, pool_size == "Pool size = 50")

levs <- round(as.numeric(as.character(new_dat$true_dynamic_prop)), 5)

new_dat$true_dynamic_prop<-as.factor(new_dat$true_dynamic_prop)
levels(new_dat$true_dynamic_prop)<-levs


#check if the cis include the true prevlance 
new_dat$within_ct_dynamic=between(new_dat$true_dynamic_p,new_dat$dynamic_q1,new_dat$dynamic_q2)*1


#get the counts 
sumrry_cts<-new_dat %>% 
  group_by(n_pools,true_dynamic_p,true_Prevalence,true_dynamic_prop) %>%
  summarise(perce_pooled=(sum(within_ct_dynamic)/n())*100,.groups = 'drop')



# new_dat$within_ct_dynamic=between(as.numeric(as.character(new_dat$true_dynamic_prop)),new_dat$prop_dynamic_q1,new_dat$prop_dynamic_q2)*1
# #get the counts 
# sumrry_cts<-new_dat %>% 
#   group_by(n_pools,true_Prevalence,true_dynamic_prop) %>%
#   summarise(perce_pooled=(sum(within_ct_dynamic)/n())*100,.groups = 'drop')



m_summary_cts=melt(sumrry_cts,id=c("n_pools","true_dynamic_p","true_Prevalence","true_dynamic_prop"))

p3<-ggplot(data=m_summary_cts,aes(x=as.factor(true_dynamic_prop),
                                  y=value,colour = as.factor(n_pools)))+
  geom_point(size=2,shape=17,alpha=1,stroke=2)+
  facet_wrap(~true_Prevalence)+
  #scale_color_viridis(discrete=T,option="D")+
  
  # scale_colour_manual(
  #  values = c("#007160", "#008e63", "#2fab63", "#78c664", "#b8e067","#f9f871")
  # )+
  xlab("True dynamic proportion")+
  ylab("Percentage (%)")+
  theme(legend.position = "bottom")+
  labs(color = "Number of pools")

p3

p1/p2/p3
#prevalence 

p2<-ggplot(data=optim_prev_esti)+
  # geom_ribbon(aes(x=n_pools,ymin=dynamic_q1,ymax=dynamic_q2))+
  geom_boxplot(aes(x=as.factor(n_pools),y=(true_static_p-optim_p_static),color=as.factor(pool_size)))+
  facet_grid(as.factor(pool_size)~as.factor(true_Prevalence))+
  # geom_abline(slope = 0, intercept = 0,        # x = y line
  #         color = "red", linetype = "dashed")+
  default_theme+
  theme(legend.position = "none")


p2

p1|p2

p2<-ggplot()+
  geom_point(data=optim_prev_esti,aes(x=true_static_p,y=optim_p_static))+
  geom_abline(slope = 1, intercept = 0,        # x = y line
              color = "red", linetype = "dashed")

p1|p2  
