library(dplyr)
library(purrr)
library(bbmle)
library(MASS)
library(tidyr)
library(future)
library(furrr)

# system.time({
# # ----- helper functions used later -----
source("~/Documents/GitHub/west_nile_virus_abm/3_figure_generation/functions/helper_function_to_test_prev.R")

#generate random observed Ct vectors 
all_data<-readRDS("synthetic_data/all_data.rds")

#calculate the prevalances: 
all_f_kde_joint=NULL
all_f_kde_joint=lapply(1:50,function(x) readRDS(file = paste0("pre_calculations/joint_kde_calc/joint_fkde_",x,".rds")))

#generate a known prevlanace:
#take pool size to be 50 
pool_sizes=50 
n_pools<-c(10,30,50,100)
# n_pools<-5

# #get the proprtion sample 
# #to have re
# p_static=c(0.001,0.01,.1)
# p_dynamic=c(0.001,0.01,.1)


p_dynamic <- c(0.001, 0.01, 0.1)
prop      <- c(0.2, 0.4, 0.6, 0.8, 0.9) #proportion for dynamic


# px_grid<-expand.grid(p_dynamic,prop)
# 
# p_static=(px_grid[,1]/px_grid[,2])-px_grid[,1]

grids <- expand.grid(p_dynamic,prop, pool_sizes, n_pools)
grids$static<-(grids$Var1/grids$Var2)-grids$Var1
#grid <- expand.grid(p_static,px_grid[,1], pool_sizes, n_pools)

grid<-data.frame(grids$static,grids$Var1,grids$Var3,grids$Var4)
# grid <- expand.grid(p_static, p_dynamic, pool_sizes, n_pools)
p_combi <- crossing(replicate = 1:100, grid)
 #p_combi<-p_combi[1:10,]
# ----- tidy up p_combi column names (use column positions because p_combi was made with crossing/expand.grid) -----
# expected column order: rep, p_static, p_dynamic, pool_size, n_pools
p_combi_tidy <- p_combi %>%
  transmute(
    repn      = replicate,
    p_static = as.numeric(.[[2]]),
    p_dynamic= as.numeric(.[[3]]),
    pool_size= as.numeric(.[[4]]),
    n_pools  = as.numeric(.[[5]])
  )

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

n_prev<-nrow(p_combi)
n_prev

ct_threshold=40


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



# helper functions & objects assumed in scope:
# source(".../helper_function_to_test_prev.R")
# ct_vec_obs(), run_one_fit(), make_pd(), is_pd_chol(), etc.
# all_f_kde_joint, all_combos, all_data, viral_loads_vec_static, viral_loads_vec_dynamic, p_combi_tidy

# ---- parallel plan ----
# choose workers; adjust if you want fewer/more cores
workers <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = workers)

# set a reproducible seed for parallel operations (change number for different RNG stream)
seed_val <- 12345
furrr_opts <- furrr::furrr_options(seed = seed_val)

# ---- (A) parallel generation of all_pools ----
# Using future_pmap so RNG inside ct_vec_obs is reproducible and parallelised
all_pools <- future_pmap(
  p_combi_tidy,
  .f = function(repn, p_static, p_dynamic, pool_size, n_pools) {
    prev_vec <- c(p_static, p_dynamic)
    ct_vec_obs(prev = prev_vec,
               pool_sizes = pool_size,
               n_pools = n_pools,
               viral_loads_vec_static = viral_loads_vec_static,
               viral_loads_vec_dynamic = viral_loads_vec_dynamic)
  },
  .options = furrr_opts
)

# ---- sanity check ----
if (length(all_pools) != nrow(p_combi_tidy)) {
  stop("all_pools length mismatch with p_combi_tidy; aborting")
}

# ---- (B) parallel MLE fits ----
nn <- nrow(p_combi_tidy)
S  <- 5000

index_list   <- seq_len(nn)
p_rows_split <- split(p_combi_tidy, seq_len(nn))  # list of one-row tibbles to pass safely
observed_list <- all_pools

# Use future_map2 to parallelise run_one_fit calls
optim_prev_esti_list <- future_map2(
  .x = index_list,
  .y = p_rows_split,
  .f = function(idx, p_row_df) {
    # run_one_fit should be defined in your environment (from earlier code)
    run_one_fit(
      index = idx,
      p_row = p_row_df,
      observed_ct_vec = observed_list[[idx]],
      S = S,
      all_f_kde_joint = all_f_kde_joint,
      all_combos = all_combos
    )
  },
  .options = furrr_opts
)

# bind results and add true probabilities (order preserved)
optim_prev_esti <- dplyr::bind_rows(optim_prev_esti_list) %>%
  dplyr::bind_cols(
    p_combi_tidy %>% dplyr::select(p_static_true = p_static, p_dynamic_true = p_dynamic)
  )

# })

saveRDS(optim_prev_esti, "optim_prev_estims_2.rds")
