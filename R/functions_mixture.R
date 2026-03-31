## functions_mixture.R

library(tidyverse)

# =========================================================
# Problem 1: Gibbs sampler
# =========================================================
gibbs_mixnorm <- function(y, K, sigma2, n_iter = 5000, burn_in = 1000,
                          init_mu = NULL, init_c = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(y)
  
  if (is.null(init_mu)) {
    probs <- seq(0.1, 0.9, length.out = K)
    mu <- as.numeric(quantile(y, probs = probs))
  } else {
    mu <- init_mu
  }
  
  if (is.null(init_c)) {
    dist_mat <- sapply(mu, function(m) abs(y - m))
    c_vec <- max.col(-dist_mat)
  } else {
    c_vec <- init_c
  }
  
  mu_draws <- matrix(NA_real_, nrow = n_iter, ncol = K)
  c_draws  <- matrix(NA_integer_, nrow = n_iter, ncol = n)
  
  relabel_by_mu <- function(mu, c_vec) {
    ord <- order(mu)
    mu_new <- mu[ord]
    new_label <- integer(length(mu))
    new_label[ord] <- seq_along(ord)
    c_new <- new_label[c_vec]
    list(mu = mu_new, c = c_new)
  }
  
  tmp <- relabel_by_mu(mu, c_vec)
  mu <- tmp$mu
  c_vec <- tmp$c
  
  for (iter in seq_len(n_iter)) {
    for (i in seq_len(n)) {
      log_prob <- dnorm(y[i], mean = mu, sd = 1, log = TRUE)
      log_prob <- log_prob - max(log_prob)
      prob <- exp(log_prob)
      prob <- prob / sum(prob)
      c_vec[i] <- sample.int(K, size = 1, prob = prob)
    }
    
    for (k in seq_len(K)) {
      idx <- which(c_vec == k)
      n_k <- length(idx)
      sum_yk <- sum(y[idx])
      
      post_var  <- 1 / (n_k + 1 / sigma2)
      post_mean <- post_var * sum_yk
      
      mu[k] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
    }
    
    tmp <- relabel_by_mu(mu, c_vec)
    mu <- tmp$mu
    c_vec <- tmp$c
    
    mu_draws[iter, ] <- mu
    c_draws[iter, ]  <- c_vec
  }
  
  keep <- (burn_in + 1):n_iter
  mu_post <- mu_draws[keep, , drop = FALSE]
  c_post  <- c_draws[keep, , drop = FALSE]
  
  mu_hat <- colMeans(mu_post)
  
  c_hat <- apply(c_post, 2, function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  })
  
  z_prob <- matrix(0, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    z_prob[, k] <- colMeans(c_post == k)
  }
  
  list(
    mu_draws = mu_draws,
    c_draws = c_draws,
    mu_hat = mu_hat,
    c_hat = c_hat,
    z_prob = z_prob,
    burn_in = burn_in
  )
}

# =========================================================
# Problem 2: Variational Bayes
# =========================================================
vb_mixnorm <- function(y, K, sigma2, max_iter = 1000, tol = 1e-8,
                       init_m = NULL, init_r = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(y)
  
  if (is.null(init_m)) {
    probs <- seq(0.1, 0.9, length.out = K)
    m <- as.numeric(quantile(y, probs = probs))
  } else {
    m <- init_m
  }
  
  s2 <- rep(1, K)
  
  if (is.null(init_r)) {
    log_r <- sapply(seq_len(K), function(k) dnorm(y, mean = m[k], sd = 1, log = TRUE))
    log_r <- sweep(log_r, 1, apply(log_r, 1, max), "-")
    r <- exp(log_r)
    r <- r / rowSums(r)
  } else {
    r <- init_r
  }
  
  elbo <- numeric(max_iter)
  
  relabel_vb <- function(m, s2, r) {
    ord <- order(m)
    list(
      m = m[ord],
      s2 = s2[ord],
      r = r[, ord, drop = FALSE]
    )
  }
  
  compute_elbo <- function(y, r, m, s2, sigma2) {
    n <- length(y)
    K <- length(m)
    
    term_mu_prior <-
      - K / 2 * log(2 * pi * sigma2) -
      (1 / (2 * sigma2)) * sum(m^2 + s2)
    
    term_c_prior <- - n * log(K)
    
    quad <- outer(y, m, FUN = function(yy, mm) yy^2 - 2 * yy * mm) +
      matrix(m^2 + s2, nrow = n, ncol = K, byrow = TRUE)
    
    term_like <-
      - n / 2 * log(2 * pi) -
      0.5 * sum(r * quad)
    
    term_qmu <- sum(0.5 * log(2 * pi * exp(1) * s2))
    
    eps <- 1e-12
    term_qc <- -sum(r * log(pmax(r, eps)))
    
    term_mu_prior + term_c_prior + term_like + term_qmu + term_qc
  }
  
  for (iter in seq_len(max_iter)) {
    log_r <- matrix(NA_real_, nrow = n, ncol = K)
    for (k in seq_len(K)) {
      log_r[, k] <- y * m[k] - 0.5 * (m[k]^2 + s2[k])
    }
    
    log_r <- sweep(log_r, 1, apply(log_r, 1, max), "-")
    r <- exp(log_r)
    r <- r / rowSums(r)
    
    Nk <- colSums(r)
    s2 <- 1 / (1 / sigma2 + Nk)
    m <- s2 * colSums(r * y)
    
    tmp <- relabel_vb(m, s2, r)
    m <- tmp$m
    s2 <- tmp$s2
    r <- tmp$r
    
    elbo[iter] <- compute_elbo(y, r, m, s2, sigma2)
    
    if (iter > 1 && abs(elbo[iter] - elbo[iter - 1]) < tol) {
      elbo <- elbo[1:iter]
      break
    }
  }
  
  c_hat <- max.col(r, ties.method = "first")
  
  list(
    m = m,
    s2 = s2,
    r = r,
    c_hat = c_hat,
    mu_hat = m,
    elbo = elbo,
    n_iter = length(elbo)
  )
}

# =========================================================
# Simulation helper functions
# =========================================================
generate_mixture_data <- function(n, mu_true, weights = NULL, sd_y = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  K <- length(mu_true)
  if (is.null(weights)) weights <- rep(1 / K, K)
  
  z <- sample(seq_len(K), size = n, replace = TRUE, prob = weights)
  y <- rnorm(n, mean = mu_true[z], sd = sd_y)
  
  list(y = y, z = z)
}

order_results <- function(mu_hat, lower, upper) {
  ord <- order(mu_hat)
  list(
    mu_hat = mu_hat[ord],
    lower = lower[ord],
    upper = upper[ord]
  )
}

simulate_one <- function(n, mu_true = c(0, 5, 10, 20), sigma2 = 1000,
                         gibbs_iter = 10000, gibbs_burn = 2000,
                         vb_max_iter = 5000, vb_tol = 1e-8,
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  K <- length(mu_true)
  
  dat <- generate_mixture_data(n = n, mu_true = mu_true, sd_y = 1)
  y <- dat$y
  
  init_mu <- mu_true + rnorm(K, 0, 0.5)
  
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

summarize_simulation <- function(sim_res) {
  bias_summary <- sim_res %>%
    pivot_longer(
      cols = c(vb_mu_hat, gibbs_mu_hat),
      names_to = "method",
      values_to = "mu_hat"
    ) %>%
    mutate(method = recode(method,
                           vb_mu_hat = "VB",
                           gibbs_mu_hat = "Gibbs"),
           error = mu_hat - true_mu) %>%
    group_by(n, component, true_mu, method) %>%
    summarise(
      bias = mean(error),
      mcse_bias = sd(error) / sqrt(n()),
      .groups = "drop"
    )
  
  coverage_summary <- sim_res %>%
    pivot_longer(
      cols = c(vb_cover, gibbs_cover),
      names_to = "method",
      values_to = "cover"
    ) %>%
    mutate(method = recode(method,
                           vb_cover = "VB",
                           gibbs_cover = "Gibbs")) %>%
    group_by(n, component, true_mu, method) %>%
    summarise(
      coverage = mean(cover),
      mcse_coverage = sqrt(coverage * (1 - coverage) / n()),
      .groups = "drop"
    )
  
  time_summary <- sim_res %>%
    group_by(n) %>%
    summarise(
      vb_time = mean(vb_time),
      gibbs_time = mean(gibbs_time),
      .groups = "drop"
    ) %>%
    pivot_longer(c(vb_time, gibbs_time),
                 names_to = "method",
                 values_to = "mean_time") %>%
    mutate(method = recode(method,
                           vb_time = "VB",
                           gibbs_time = "Gibbs"))
  
  list(
    bias = bias_summary,
    coverage = coverage_summary,
    time = time_summary
  )
}

plot_bias <- function(bias_summary) {
  ggplot(bias_summary,
         aes(x = factor(component), y = bias, color = method, group = method)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(position = position_dodge(width = 0.35), size = 2) +
    geom_errorbar(
      aes(ymin = bias - 1.96 * mcse_bias,
          ymax = bias + 1.96 * mcse_bias),
      width = 0.15,
      position = position_dodge(width = 0.35)
    ) +
    facet_wrap(~ n, scales = "free_y") +
    labs(
      x = "Component",
      y = "Bias of mu_hat",
      title = "Bias of estimated component means with MCSE error bars"
    ) +
    theme_bw()
}

plot_coverage <- function(coverage_summary) {
  ggplot(coverage_summary,
         aes(x = factor(component), y = coverage, color = method, group = method)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_point(position = position_dodge(width = 0.35), size = 2) +
    geom_errorbar(
      aes(ymin = coverage - 1.96 * mcse_coverage,
          ymax = coverage + 1.96 * mcse_coverage),
      width = 0.15,
      position = position_dodge(width = 0.35)
    ) +
    facet_wrap(~ n) +
    labs(
      x = "Component",
      y = "Coverage",
      title = "Coverage of estimated component means with MCSE error bars"
    ) +
    theme_bw()
}