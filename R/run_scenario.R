## run_scenario.R

args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[1])

library(here)
source(here("R", "functions_mixture.R"))

dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)

sim_res <- run_simulation_scenario(
  n = n,
  nsim = 500,
  mu_true = c(0, 5, 10, 20),
  sigma2 = 1000,
  gibbs_iter = 10000,
  gibbs_burn = 2000,
  vb_max_iter = 5000,
  vb_tol = 1e-8,
  seed = 2026 + n
)

saveRDS(sim_res, file = here("results", paste0("sim_n", n, ".rds")))