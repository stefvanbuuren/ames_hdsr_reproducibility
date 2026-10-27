# Goals 1-3 production run: original methods plus FCS_train_fill (prediction only).
# Goals 2-3 prediction models use log_price ~ X + Z; Goal 4 (separate script) keeps ~ X.
source('ames_helpers.R')
d <- ames_preprocess('Gr_Liv_Area')
stopifnot(nrow(d) == 2925L, max(d$Gr_Liv_Area) <= 4000)

# Production configuration for Table 1. For a quick smoke test, lower nrep, m_train and m_fill.
set.seed(20261022L)
t0 <- Sys.time()
results_123 <- run_ames_comparison(
  targets = 'Gr_Liv_Area',
  nrep = 100L,
  m = 100L,
  fcs = TRUE,
  m_train = 100L,
  m_fill = 1000L
)
print(Sys.time() - t0)
summary_123 <- summary_table(results_123)
print(summary_123)
write.csv(summary_123, 'ames_goals123_FCS_production_summary.csv', row.names = FALSE)
saveRDS(results_123, 'ames_goals123_FCS_production_results.rds')

# Reproducibility metadata
writeLines(capture.output(sessionInfo()), 'sessionInfo_goals123.txt')
