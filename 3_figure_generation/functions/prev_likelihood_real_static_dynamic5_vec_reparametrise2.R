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
    
    obs_mat<-all_combos[[x]]
    
    # compute multinomial densities vectorised over rows.
    log_dens_mult<-apply(obs_mat,1,function(x) dmultinom(x, prob = prob,log = TRUE))
    
    #load the joint density 
    dens_joint<-unlist(lapply(f_kde_joint[[x]], function(f) f(yi)))
    
    # combined probability for this pool (sum over all possible (k,ks) cases)
    #p_i <- sum(dens_mult * dens_s * dens_d)
    log_sums<-log_dens_mult+log(dens_joint)
    log_p_i<-log_sum_exp(log_sums)
    
    # guard against zero (log(0) -> -Inf). Use eps floor for p_i too.
    if (!is.finite(log_p_i)) {
      log_p_i <- log(.Machine$double.eps)
    }
    

    loglik <- loglik + log_p_i
  }
  
  if (neg_lik) -loglik else loglik
}
