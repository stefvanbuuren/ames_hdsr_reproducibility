# Parallel equivalent of run_goals123.R: each replication is seeded by seed + rep inside ames_one(),
# so results are order-independent and identical to the sequential loop.
source('ames_helpers.R')
library(parallel)
t0 <- Sys.time()
res <- mclapply(1:100, function(b) ames_one('Gr_Liv_Area', b, 100L, 1/3, .25, 20261022L,
                                            'linear', fcs = TRUE, m_train = 100L, m_fill = 1000L),
                mc.cores = 12, mc.preschedule = FALSE)
stopifnot(!any(vapply(res, inherits, logical(1), 'try-error')))
results_123 <- do.call(rbind, res); rownames(results_123) <- NULL
print(Sys.time() - t0)
s <- summary_table(results_123)
write.csv(s, 'ames_goals123_FCS_production_summary.csv', row.names = FALSE)
saveRDS(results_123, 'ames_goals123_FCS_production_results.rds')
