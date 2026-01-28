prev_likelihood_real_static_dynamic5_vec_reparametrise2<- function(
    f_kde_joint,
    all_combos,
    a,
    b,
    n_pools,
    pool_sizes,
    observed_ct_vec,
    neg_lik = TRUE
) {
  #stick-breaking to ensure the sum is <1 
  p_static <- plogis(a)                       # (0,1)
  p_dynamic <- plogis(b) * (1 - p_static)     # ensures p_static + p_dynamic <= 1
  
  #x is a vector: 
  log_sum_exp <- function(x) {
    if (all(!is.finite(x))) return(-Inf)
    m <- max(x)
    m + log(sum(exp(x - m)))
  }
  

  prob <- c(p_static, p_dynamic, 1 - p_static - p_dynamic)
  eps  <- 1e-5
  loglik <- 0
  
  for (i in seq_len(n_pools)) {
    x <- pool_sizes[i]
    yi <- observed_ct_vec[i]
    
    # # build vectors of all (k, ks) combinations:
    # # for k = 0..x, ks runs 0..k (so number of rows = sum_{k=0}^x (k+1) = (x+1)(x+2)/2)
    # k_vals <- 0:x
    # reps   <- k_vals + 1L               # number of ks values for each k
    # # ks_vec: 0, 0:1, 0:2, ..., 0:x  (unlisted)
    # ks_vec <- unlist(lapply(k_vals, function(k) 0:k), use.names = FALSE)
    # # kd_vec: for each k, kd = k - ks, so gives k:(0)
    # kd_vec <- unlist(lapply(k_vals, function(k) k - (0:k)), use.names = FALSE)
    # # kn_vec: kn = x - k, repeated (k+1) times for each k
    # kn_vec <- rep(x - k_vals, times = reps)
    # 
    # 
    # # assemble observation rows for multinomial density
    # obs_mat <- cbind(ks_vec, kd_vec, kn_vec)
    obs_mat<-all_combos[[x]]
    
    # compute multinomial densities vectorised over rows.
    # dens_mult <- dmultinom(obs_mat, prob = prob)
    log_dens_mult<-apply(obs_mat,1,function(x) dmultinom(x, prob = prob,log = TRUE))
    
    
    # evaluate KDE densities for static (indexed by ks) and dynamic (indexed by kd)
    # f_kde_static[[i]] is a list of length (x+1) of functions (indexed by count+1)
    # we evaluate each appropriate function at yi and apply floor eps
    #######
    # joint_k<-joint_kde(B=100,obs_mat,viral_loads_vec_static,viral_loads_vec_dynamic)
    #
    
    dens_joint<-unlist(lapply(f_kde_joint[[x]], function(f) f(yi)))
    
    # combined probability for this pool (sum over all possible (k,ks) cases)
    #p_i <- sum(dens_mult * dens_s * dens_d)
    log_sums<-log_dens_mult+log(dens_joint)
    log_p_i<-log_sum_exp(log_sums)
    
    # guard against zero (log(0) -> -Inf). Use eps floor for p_i too.
    if (!is.finite(log_p_i)) {
      log_p_i <- log(.Machine$double.eps)
    }
    
    #loglik <- loglik + log(p_i)
    loglik <- loglik + log_p_i
  }
  
  if (neg_lik) -loglik else loglik
}
