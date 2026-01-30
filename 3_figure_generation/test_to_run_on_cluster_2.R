library(dplyr)
library(purrr)
library(bbmle)
library(MASS)
library(tidyr)


# ----- helper functions used later -----
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
#n_pools<-5

#get the proprtion sample 
#to have re
p_static=c(0.001,0.01,.1)
p_dynamic=c(0.001,0.01,.1)


grid <- expand.grid(p_static, p_dynamic, pool_sizes, n_pools)
p_combi <- crossing(rep = 1:100, grid)
# ----- tidy up p_combi column names (use column positions because p_combi was made with crossing/expand.grid) -----
# expected column order: rep, p_static, p_dynamic, pool_size, n_pools
p_combi_tidy <- p_combi %>%
  transmute(
    rep      = .[[1]],
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




# ----- 1) replace the loop that builds `all_pools` with pmap -----
# Result: list of length nrow(p_combi_tidy), each element = vector of pooled Cts for that scenario
all_pools <- p_combi_tidy %>%
  # pmap will pass columns into the anonymous function in order (rep, p_static, p_dynamic, pool_size, n_pools)
  pmap(function(rep, p_static, p_dynamic, pool_size, n_pools) {
    prev_vec <- c(p_static, p_dynamic)
    ct_vec_obs(prev = prev_vec,
               pool_sizes = pool_size,
               n_pools = n_pools,
               viral_loads_vec_static = viral_loads_vec_static,
               viral_loads_vec_dynamic = viral_loads_vec_dynamic)
  })

# all_pools is now a list; you used it the same way later (all_pools[[i]])




# ----- 2) replace the final for-loop (MLE over each prevalence row) with map2 / pmap -----
nn <- nrow(p_combi_tidy)

S <- 5000
# Build a list of indices and observed_ct_vec (to feed safely into mapping function)
index_list <- seq_len(nn)
observed_list <- all_pools  # list of CT vectors, same length as p_combi_tidy

# run fits in series using purrr::map2 (index, row) -> returns list of tibbles, then bind_rows
optim_prev_esti_list <- purrr::map2(
  .x = index_list,
  .y = split(p_combi_tidy, seq_len(nn)), # each element is a one-row data.frame
  .f = function(idx, p_row_df) {
    run_one_fit(
      index = idx,
      p_row = p_row_df,
      observed_ct_vec = observed_list[[idx]],
      S = S,
      all_f_kde_joint = all_f_kde_joint,
      all_combos = all_combos
    )
  }
)

optim_prev_esti <- dplyr::bind_rows(optim_prev_esti_list)

optim_prev_esti <- dplyr::bind_rows(optim_prev_esti_list) %>%
  dplyr::bind_cols(
    p_combi_tidy %>% dplyr::select(p_static_true = p_static, p_dynamic_true = p_dynamic)
  )


# ---- optim_prev_esti is the same shape as earlier rbind results ----
