# Table 2 (Goal 4): repeated sampling from the empirical Ames population.
# Source files from the project root. Outputs include the PMM sensitivity analysis.
source('ames_helpers.R')
source('ames_goal4.R')

# Production configuration: 1,000 bootstrap-like samples, 1,000 houses each.
# For a smoke test, temporarily set nrep = 20L.
results_goal4 <- run_ames_goal4(
  nrep = 1000L,
  n_sample = 1000L,
  m = 20L,
  seed = 20261022L
)
summary_goal4 <- summary_ames_goal4(results_goal4)
print(summary_goal4)
write.csv(summary_goal4, 'ames_goal4_summary.csv', row.names = FALSE)
saveRDS(results_goal4, 'ames_goal4_results.rds')
writeLines(capture.output(sessionInfo()), 'sessionInfo_goal4.txt')
