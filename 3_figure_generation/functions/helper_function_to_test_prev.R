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


make_pd <- function(Sigma) {
  S <- (Sigma + t(Sigma)) / 2
  if (is_pd_chol(S)) return(S)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for nearPD(). Please install it.")
  }
  as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
}


is_pd_chol <- function(Sigma) {
  out <- try(chol(Sigma), silent = TRUE)
  !inherits(out, "try-error")
}

# wrapper to run single fit and bootstrap -> returns a data.frame or NA-row on failure
run_one_fit <- function(index, p_row, observed_ct_vec, S = 5000, all_f_kde_joint, all_combos) {
  # index: integer (for week labeling)
  # p_row: one-row tibble / list with pool_size and n_pools etc.
  tryCatch({
    b_pools <- as.integer(p_row$n_pools)
    sample_size_i <- rep(as.integer(p_row$pool_size), b_pools)
    
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
    a_sim <- boots[, 1]
    b_sim <- boots[, 2]
    
    p_static_sim  <- plogis(a_sim)
    p_dynamic_sim <- plogis(b_sim) * (1 - p_static_sim)
    
    p_dynamic_sim_prop <- p_dynamic_sim / (p_dynamic_sim + p_static_sim)
    p_static_sim_prop  <- p_static_sim / (p_dynamic_sim + p_static_sim)
    
    ci_p_static  <- quantile(p_static_sim,  probs = c(.025, .975), na.rm = TRUE)
    ci_p_dynamic <- quantile(p_dynamic_sim, probs = c(.025, .975), na.rm = TRUE)
    
    ci_prop_static  <- quantile(p_static_sim_prop,  probs = c(.025, .975), na.rm = TRUE)
    ci_prop_dynamic <- quantile(p_dynamic_sim_prop,  probs = c(.025, .975), na.rm = TRUE)
    
    tibble::tibble(
      week            = index,
      optim_p_static  = p_static_hat,
      optim_p_dynamic = p_dynamic_hat,
      static_q1       = ci_p_static[1],
      static_q2       = ci_p_static[2],
      dynamic_q1      = ci_p_dynamic[1],
      dynamic_q2      = ci_p_dynamic[2],
      prop_dynamic_q1 = ci_prop_dynamic[1],
      prop_dynamic_q2 = ci_prop_dynamic[2],
      prop_static_q1  = ci_prop_static[1],
      prop_static_q2  = ci_prop_static[2]
    )
  }, error = function(e) {
    # on error return a row of NAs but keep the week index so result length stays constant
    warning(sprintf("fit failed for index %s: %s", index, e$message))
    tibble::tibble(
      week            = index,
      optim_p_static  = NA_real_,
      optim_p_dynamic = NA_real_,
      static_q1       = NA_real_,
      static_q2       = NA_real_,
      dynamic_q1      = NA_real_,
      dynamic_q2      = NA_real_,
      prop_dynamic_q1 = NA_real_,
      prop_dynamic_q2 = NA_real_,
      prop_static_q1  = NA_real_,
      prop_static_q2  = NA_real_
    )
  })
}

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

