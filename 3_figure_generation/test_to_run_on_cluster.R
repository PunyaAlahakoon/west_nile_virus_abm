

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



# parallel version using foreach + doParallel + doRNG (works on Windows & *nix)
library(doParallel)
library(foreach)
library(doRNG)   # ensures reproducible RNG across workers
library(bbmle)   # required by your mle2 call
library(MASS)    # for mvrnorm

# number of workers: adjust as desired (detectCores()-1 is common)
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# set a reproducible seed for the whole job
registerDoRNG(seed = 12345)

# Ensure these objects/functions are available on the workers:
# - prev_likelihood_real_static_dynamic5_vec_reparametrise2
# - make_pd
# - all_f_kde_joint, all_combos, all_pools, p_combi, S, nn
# If they are in your global env, export them (add names if different)
exports <- c("prev_likelihood_real_static_dynamic5_vec_reparametrise2",
             "make_pd",
             "all_f_kde_joint","all_combos","all_pools","p_combi","S","nn")

# Run parallel foreach
optim_prev_esti <- foreach(i = 1:2,
                           .combine = rbind,
                           .packages = c("bbmle","MASS"),   # packages needed inside
                           .export = exports) %dopar% {
                             
                          
                        
                               b_pools <- as.numeric(p_combi[i,5])
                               sample_size_i <- as.numeric(rep(p_combi[i,4], b_pools))
                               observed_ct_vec <- all_pools[[i]]
                               
                               fit_a <- bbmle::mle2(
                                 prev_likelihood_real_static_dynamic5_vec_reparametrise2,
                                 start = list(a = qlogis(0.005), b = qlogis(0.005)),
                                 data = list(
                                   f_kde_joint = all_f_kde_joint,
                                   all_combos = all_combos,
                                   n_pools = b_pools,
                                   pool_sizes = sample_size_i,
                                   observed_ct_vec = observed_ct_vec,
                                   neg_lik = TRUE
                                 ),
                                 method = "L-BFGS-B",
                                 lower = c(-50, -50),
                                 upper = c(50, 50)
                               )
                               
                               a_hat <- coef(fit_a)[["a"]]
                               b_hat <- coef(fit_a)[["b"]]
                               p_static_hat <- plogis(a_hat)
                               p_dynamic_hat <- plogis(b_hat) * (1 - p_static_hat)
                               
                               cov_mat <- make_pd(vcov(fit_a))
                               boots <- MASS::mvrnorm(S, mu = coef(fit_a), Sigma = cov_mat)
                               a_sim <- boots[,1]
                               b_sim <- boots[,2]
                               
                               p_static_sim  <- plogis(a_sim)
                               p_dynamic_sim <- plogis(b_sim) * (1 - p_static_sim)
                               
                               p_dynamic_sim_prop <- p_dynamic_sim / (p_dynamic_sim + p_static_sim)
                               p_static_sim_prop  <- p_static_sim  / (p_dynamic_sim + p_static_sim)
                               
                               ci_p_static  <- quantile(p_static_sim,  probs = c(.025, .975))
                               ci_p_dynamic <- quantile(p_dynamic_sim, probs = c(.025, .975))
                               
                               ci_prop_static  <- quantile(p_static_sim_prop,  probs = c(.025, .975), na.rm = TRUE)
                               ci_prop_dynamic <- quantile(p_dynamic_sim_prop,  probs = c(.025, .975), na.rm = TRUE)
                               
                               # Return a single-row data.frame 
                               data.frame(
                                 week = i,
                                 optim_p_static   = p_static_hat,
                                 optim_p_dynamic  = p_dynamic_hat,
                                 static_q1        = ci_p_static[1],
                                 static_q2        = ci_p_static[2],
                                 dynamic_q1       = ci_p_dynamic[1],
                                 dynamic_q2       = ci_p_dynamic[2],
                                 prop_dynamic_q1  = ci_prop_dynamic[1],
                                 prop_dynamic_q2  = ci_prop_dynamic[2],
                                 prop_static_q1   = ci_prop_static[1],
                                 prop_static_q2   = ci_prop_static[2]
                               )
                               
                      
                           }



# stop cluster when done
stopCluster(cl)
registerDoSEQ()  # return to sequential backend if needed



for (i in 1:nn) {
  b_pools<-as.numeric(p_combi[i,5])
  sample_size_i<-as.numeric(rep(p_combi[i,4],b_pools))
  #observed_ct_vec<-as.numeric(all_pools[i,])
  observed_ct_vec<-all_pools[[i]]
 
  
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

    )
  )
  
  
}


saveRDS(optim_prev_esti, "optim_prev_esti_5.rds")
