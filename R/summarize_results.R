## summarize_results.R
library(tidyverse)
library(here)
source(here("R", "functions_mixture.R"))

res100 <- readRDS(here("results", "sim_n100.rds"))
res1000 <- readRDS(here("results", "sim_n1000.rds"))
res10000 <- readRDS(here("results", "sim_n10000.rds"))

all_res <- bind_rows(res100, res1000, res10000)

summ <- summarize_simulation(all_res)

bias_summary <- summ$bias
coverage_summary <- summ$coverage
time_summary <- all_res %>%
  group_by(n) %>%
  summarise(
    mean_vb_time = mean(vb_time),
    mean_gibbs_time = mean(gibbs_time),
    .groups = "drop"
  )

print(bias_summary)
print(coverage_summary)
print(time_summary)

p1 <- plot_bias(bias_summary)
p2 <- plot_coverage(coverage_summary)

ggsave(here("plots", "bias_plot.png"), p1, width = 10, height = 6)
ggsave(here("plots", "coverage_plot.png"), p2, width = 10, height = 6)