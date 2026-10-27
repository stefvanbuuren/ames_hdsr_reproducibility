# Ames Goal 4: repeated sampling and confidence-interval coverage
# Version 2026-10-22
# Uses 2925 Ames houses with Gr_Liv_Area <= 4000 as an EMPIRICAL population.
# Population estimand: least-squares linear-projection slope of log(Sale_Price)
# on Gr_Liv_Area, with an intercept. This is not a causal effect.
# Draw independent samples WITH replacement, mask 1/3 of Gr_Liv_Area MCAR,
# and evaluate bias, sampling SD, mean reported SE, 95% CI coverage and width.
# Imputation specifications mirror ames_helpers.R:
#   Z = log_lot_area + Lot_Config + Bldg_Type
#   LM_mean_Z, LM_mean_ZY, LM_draw_Z, LM_draw_ZY, LM_MI_ZY, LM_MI_PMM_ZY
# CompleteData provides a benchmark; it is not an imputation method.
# All within-analysis variances use heteroskedasticity-robust HC3.
# For MI, combine HC3 variances with Rubin's rules (mice::pool.scalar).
# Caveat: treating the 2925 retained houses as a population means results are conditional
# on this empirical distribution, not a population-level housing claim.
# Requires: modeldata, mice, sandwich, and the accompanying coverage-by-status script.
# Uses the shared helpers in ames_helpers.R.

ames_goal4_load <- function(helper = 'ames_helpers.R') {
  if (!exists('ames_preprocess', mode='function') ||
      !exists('ames_mi_train', mode='function')) {
    if (!file.exists(helper)) stop('Place ',helper,' in your working directory, or source it first.')
    source(helper, local=.GlobalEnv)
  }
  invisible(TRUE)
}

ames_goal4_ci <- function(fit, term='Gr_Liv_Area', level=.95) {
  if (!requireNamespace('sandwich',quietly=TRUE)) stop('Install package sandwich')
  estimate <- unname(stats::coef(fit)[term])
  if (length(estimate)!=1L || !is.finite(estimate)) stop('Coefficient not found: ',term)
  se <- sqrt(unname(sandwich::vcovHC(fit,type='HC3')[term,term]))
  # Large-sample robust Wald interval. Identical convention across methods.
  critical <- stats::qnorm((1+level)/2)
  c(estimate=estimate, se=se, lower=estimate-critical*se,
    upper=estimate+critical*se)
}

ames_goal4_mi_ci <- function(completed, formula, term='Gr_Liv_Area', level=.95) {
  if (!requireNamespace('mice',quietly=TRUE)) stop('Install package mice')
  if (!requireNamespace('sandwich',quietly=TRUE)) stop('Install package sandwich')
  fits <- lapply(completed,function(d) stats::lm(formula,data=d))
  q <- vapply(fits,function(fit) unname(stats::coef(fit)[term]),numeric(1))
  u <- vapply(fits,function(fit)
    unname(sandwich::vcovHC(fit,type='HC3')[term,term]),numeric(1))
  if (any(!is.finite(q)) || any(!is.finite(u)) || any(u<0))
    stop('Non-finite estimate or HC3 variance in MI')
  # Rubin's rules with HC3 as the within-imputation variance.
  # n and k supply the complete-data degrees of freedom for Barnard-Rubin.
  n <- nrow(completed[[1]])
  pooled <- mice::pool.scalar(Q=q,U=u,n=n,k=length(stats::coef(fits[[1]])))
  se <- sqrt(pooled$t)
  critical <- stats::qt((1+level)/2,df=pooled$df)
  c(estimate=pooled$qbar,se=se,
    lower=pooled$qbar-critical*se,upper=pooled$qbar+critical*se)
}

# Predictive mean matching sensitivity analysis: impute X conditional on Z and Y.
# Only Gr_Liv_Area is imputed; all other variables are retained unchanged.
ames_mi_pmm <- function(incomplete, target, z, m=20L, seed=123L) {
  if (!requireNamespace('mice',quietly=TRUE)) stop('Install package mice')
  vars <- unique(c(target,z,'log_price'))
  dat <- incomplete[,vars,drop=FALSE]
  pred <- mice::make.predictorMatrix(dat)
  pred[,] <- 0
  pred[target,setdiff(vars,target)] <- 1
  meth <- mice::make.method(dat)
  meth[] <- ''
  meth[target] <- 'pmm'
  imp <- mice::mice(dat,m=m,method=meth,predictorMatrix=pred,
                    donors=5,maxit=5,seed=seed,printFlag=FALSE)
  lapply(seq_len(m),function(i) {
    completed <- incomplete
    completed[[target]] <- mice::complete(imp,i)[[target]]
    completed
  })
}

ames_goal4_one <- function(population, population_slope, rep,
                           n_sample=1000L, missing_fraction=1/3,
                           m=20L, seed=20261022L, level=.95) {
  target <- 'Gr_Liv_Area'
  z <- attr(population,'z')
  if (is.null(z)) z <- c('log_lot_area','Lot_Config','Bldg_Type')
  # Replication-specific, reproducible sampling/masking.
  set.seed(seed + rep)
  ids <- sample.int(nrow(population), size=n_sample, replace=TRUE)
  truth <- population[ids,,drop=FALSE]
  rownames(truth) <- NULL
  mask <- sample.int(n_sample,size=round(n_sample*missing_fraction),replace=FALSE)
  incomplete <- truth
  incomplete[[target]][mask] <- NA_real_

  fit_z <- ames_fit(incomplete,target,z,use_y=FALSE)
  fit_zy <- ames_fit(incomplete,target,z,use_y=TRUE)
  mean_z <- ames_fill(incomplete,fit_z,target,stochastic=FALSE)
  mean_zy <- ames_fill(incomplete,fit_zy,target,stochastic=FALSE)
  set.seed(seed + 10000L*rep + 1L)
  draw_z <- ames_fill(incomplete,fit_z,target,stochastic=TRUE)
  set.seed(seed + 10000L*rep + 2L)
  draw_zy <- ames_fill(incomplete,fit_zy,target,stochastic=TRUE)
  mi <- ames_mi_train(incomplete,target,z,m=m,seed=seed + 10000L*rep + 3L)
  mi_pmm <- ames_mi_pmm(incomplete,target,z,m=m,seed=seed + 10000L*rep + 4L)

  f <- stats::reformulate(target,response='log_price')
  methods <- list(CompleteData=truth, LM_mean_Z=mean_z,
                  LM_mean_ZY=mean_zy, LM_draw_Z=draw_z,
                  LM_draw_ZY=draw_zy)
  out <- lapply(names(methods),function(method) {
    ci <- ames_goal4_ci(stats::lm(f,data=methods[[method]]),target,level)
    data.frame(rep=rep,method=method,estimate=ci['estimate'],
               reported_se=ci['se'],ci_lower=ci['lower'],ci_upper=ci['upper'],
               stringsAsFactors=FALSE,row.names=NULL)
  })
  ci <- ames_goal4_mi_ci(mi,f,target,level)
  out[[length(out)+1L]] <- data.frame(
    rep=rep,method='LM_MI_ZY',estimate=ci['estimate'],
    reported_se=ci['se'],ci_lower=ci['lower'],ci_upper=ci['upper'],
    stringsAsFactors=FALSE,row.names=NULL)
  ci <- ames_goal4_mi_ci(mi_pmm,f,target,level)
  out[[length(out)+1L]] <- data.frame(
    rep=rep,method='LM_MI_PMM_ZY',estimate=ci['estimate'],
    reported_se=ci['se'],ci_lower=ci['lower'],ci_upper=ci['upper'],
    stringsAsFactors=FALSE,row.names=NULL)
  ans <- do.call(rbind,out)
  ans$population_slope <- population_slope
  ans$covered <- as.integer(ans$ci_lower <= population_slope &
                            ans$ci_upper >= population_slope)
  ans$ci_width <- ans$ci_upper-ans$ci_lower
  ans$n_sample <- n_sample
  ans$n_masked <- length(mask)
  ans
}

run_ames_goal4 <- function(nrep=100L,n_sample=1000L,m=20L,
                           missing_fraction=1/3,seed=20261022L,
                           level=.95,verbose=TRUE,
                           helper='ames_helpers.R') {
  ames_goal4_load(helper)
  stopifnot(nrep>=2L,n_sample>=10L,m>=2L,
            missing_fraction>0,missing_fraction<1,
            level>0,level<1)
  pop <- ames_preprocess('Gr_Liv_Area')
  pop <- pop[pop$Gr_Liv_Area <= 4000, , drop=FALSE]
  if (nrow(pop)!=2925L) stop('Expected 2925 houses after restriction; got ',nrow(pop))
  beta <- ames_slope(pop,'Gr_Liv_Area')
  if (verbose) message('Empirical population N=',nrow(pop),
                       '; slope=',signif(beta,7))
  res <- vector('list',nrep)
  for (b in seq_len(nrep)) {
    if (verbose && (b==1L || b%%10L==0L || b==nrep))
      message('Goal 4 replication ',b,'/',nrep)
    res[[b]] <- ames_goal4_one(pop,beta,b,n_sample,
                               missing_fraction,m,seed,level)
  }
  ans <- do.call(rbind,res)
  rownames(ans) <- NULL
  attr(ans,'population_slope') <- beta
  ans
}

summary_ames_goal4 <- function(results) {
  split_rows <- split(results,results$method)
  out <- lapply(split_rows,function(d) {
    theta <- d$population_slope[1]
    data.frame(method=d$method[1],
               nrep=nrow(d),
               population_slope=theta,
               mean_estimate=mean(d$estimate),
               bias=mean(d$estimate-theta),
               relative_bias_pct=100*mean(d$estimate-theta)/theta,
               empirical_sd=stats::sd(d$estimate),
               mean_reported_se=mean(d$reported_se),
               se_ratio=mean(d$reported_se)/stats::sd(d$estimate),
               ci_coverage=mean(d$covered),
               ci_coverage_mcse=sqrt(mean(d$covered)*(1-mean(d$covered))/nrow(d)),
               mean_ci_width=mean(d$ci_width),row.names=NULL)
  })
  ans <- do.call(rbind,out)
  ans <- ans[match(c('CompleteData','LM_mean_Z','LM_mean_ZY',
                     'LM_draw_Z','LM_draw_ZY','LM_MI_ZY','LM_MI_PMM_ZY'),ans$method),]
  rownames(ans) <- NULL
  ans
}

# Usage:
# source('ames_helpers.R')
# source('ames_goal4_hc3_restricted.R')
# pilot <- run_ames_goal4(nrep=20, n_sample=1000, m=20)
# summary_ames_goal4(pilot)
# results <- run_ames_goal4(nrep=1000, n_sample=1000, m=20)
# summary_ames_goal4(results)
# write.csv(results,'ames_goal4_replicates.csv',row.names=FALSE)
# write.csv(summary_ames_goal4(results),'ames_goal4_summary.csv',row.names=FALSE)
