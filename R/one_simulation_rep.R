simulate_one <- function(n, mu_true = c(0, 5, 10, 20), sigma2 = 1000,
                         gibbs_iter = 10000, gibbs_burn = 2000,
                         vb_max_iter = 5000, vb_tol = 1e-8,
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  K <- length(mu_true)
  
  # -----------------------------
  # Generate data
  # -----------------------------
  dat <- generate_mixture_data(n = n, mu_true = mu_true, sd_y = 1)
  y <- dat$y
  
  # sensible initialization near true values
  init_mu <- mu_true + rnorm(K, 0, 0.5)
  
  # -----------------------------
  # Variational Bayes
  # -----------------------------
  time_vb <- system.time({
    fit_vb <- vb_mixnorm(
      y = y,
      K = K,
      sigma2 = sigma2,
      max_iter = vb_max_iter,
      tol = vb_tol,
      init_m = init_mu
    )
  })["elapsed"]
  
  vb_mu_hat <- fit_vb$mu_hat
  vb_lower <- fit_vb$m - 1.96 * sqrt(fit_vb$s2)
  vb_upper <- fit_vb$m + 1.96 * sqrt(fit_vb$s2)
  
  vb_ord <- order_results(vb_mu_hat, vb_lower, vb_upper)
  
  # -----------------------------
  # Gibbs sampler
  # -----------------------------
  time_gibbs <- system.time({
    fit_gibbs <- gibbs_mixnorm(
      y = y,
      K = K,
      sigma2 = sigma2,
      n_iter = gibbs_iter,
      burn_in = gibbs_burn,
      init_mu = init_mu
    )
  })["elapsed"]
  
  post_draws <- fit_gibbs$mu_draws[(gibbs_burn + 1):gibbs_iter, , drop = FALSE]
  
  gibbs_mu_hat <- colMeans(post_draws)
  gibbs_lower <- apply(post_draws, 2, quantile, probs = 0.025)
  gibbs_upper <- apply(post_draws, 2, quantile, probs = 0.975)
  
  gibbs_ord <- order_results(gibbs_mu_hat, gibbs_lower, gibbs_upper)
  
  tibble(
    n = n,
    component = seq_len(K),
    
    true_mu = mu_true,
    
    vb_mu_hat = vb_ord$mu_hat,
    vb_lower = vb_ord$lower,
    vb_upper = vb_ord$upper,
    vb_cover = as.integer(mu_true >= vb_ord$lower & mu_true <= vb_ord$upper),
    vb_time = as.numeric(time_vb),
    
    gibbs_mu_hat = gibbs_ord$mu_hat,
    gibbs_lower = gibbs_ord$lower,
    gibbs_upper = gibbs_ord$upper,
    gibbs_cover = as.integer(mu_true >= gibbs_ord$lower & mu_true <= gibbs_ord$upper),
    gibbs_time = as.numeric(time_gibbs)
  )
}

run_simulation_scenario <- function(n, nsim = 500,
                                    mu_true = c(0, 5, 10, 20),
                                    sigma2 = 1000,
                                    gibbs_iter = 10000,
                                    gibbs_burn = 2000,
                                    vb_max_iter = 5000,
                                    vb_tol = 1e-8,
                                    seed = 2026) {
  
  res_list <- vector("list", nsim)
  
  for (s in seq_len(nsim)) {
    if (s %% 25 == 0) message("n = ", n, ", sim = ", s, "/", nsim)
    
    res_list[[s]] <- simulate_one(
      n = n,
      mu_true = mu_true,
      sigma2 = sigma2,
      gibbs_iter = gibbs_iter,
      gibbs_burn = gibbs_burn,
      vb_max_iter = vb_max_iter,
      vb_tol = vb_tol,
      seed = seed + s
    ) %>%
      mutate(sim = s)
  }
  
  bind_rows(res_list)
}
