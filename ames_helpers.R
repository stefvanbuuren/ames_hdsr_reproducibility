# Ames Housing study: shared helpers for Goals 1-4.
# Source this file from run_goals123.R or run_goal4.R.
# Data: modeldata::ames; Gr_Liv_Area <= 4000 (N = 2925).
# Dependencies: modeldata, mice (development version with train/fill for Goal 3 FCS).
# Original audited simulation implementation retained; only documentation and
# entry-point naming have been cleaned for public release.
#
# 1/3 MCAR masking applies to cases with observed target, independently in train/test.
# All goal-2/3 test outcomes are strictly held out until scoring.
# FCS_train_fill uses stored X|Y,Z and Y|X,Z models for prediction only.
# Goal 3 uses a mixture of predictive normal distributions, NOT averaged PI endpoints.
# These are plug-in predictive distributions: regression parameter uncertainty and
# uncertainty in estimated X|Z model parameters are NOT fully integrated.
# Near-nominal empirical coverage does NOT establish theoretical calibration.
# Oracle provides a benchmark using the true X even for artificially masked cases.
# Goals 2-3 (prediction) use log_price ~ X + Z for every method, matching
# FCS_train_fill, which predicts Y from X and Z. Goal 4 keeps the inferential
# estimand unchanged: the slope of log_price ~ X (see ames_slope()).
# Goal 4 slope discrepancy is relative to a complete-data fitted slope, NOT bias.
# Replications resplit and remask a fixed Ames dataset; not population sampling.

AMES_CANDIDATES <- c('Lot_Frontage', 'Gr_Liv_Area', 'Total_Bsmt_SF')
AMES_METHODS <- c(
  'LM_mean_Z',
  'LM_mean_ZY',
  'LM_draw_Z',
  'LM_draw_ZY',
  'LM_MI_ZY',
  'Oracle',
  'FCS_train_fill'
)

ames_preprocess <- function(target) {
  stopifnot(target %in% AMES_CANDIDATES)
  if (!requireNamespace('modeldata', quietly = TRUE)) {
    stop('Install modeldata')
  }
  a <- as.data.frame(modeldata::ames)
  if (
    !all(
      c(
        'Sale_Price',
        'Lot_Area',
        'Lot_Config',
        'Bldg_Type',
        'Gr_Liv_Area',
        'Total_Bsmt_SF',
        'Lot_Frontage'
      ) %in%
        names(a)
    )
  ) {
    stop('Unexpected modeldata::ames column names')
  }
  a$Lot_Frontage <- as.numeric(a$Lot_Frontage)
  a$Lot_Frontage[a$Lot_Frontage == 0] <- NA_real_ # Ames sentinel only
  a$Gr_Liv_Area <- as.numeric(a$Gr_Liv_Area)
  a$Total_Bsmt_SF <- as.numeric(a$Total_Bsmt_SF) # legitimate zeros retained
  a$log_price <- log(a$Sale_Price)
  a$log_lot_area <- log(a$Lot_Area)
  # Shared, pre-specified auxiliary predictors across all three candidate targets.
  # For living/basement area, consider alternative Z sets as a sensitivity analysis.
  z <- c('log_lot_area', 'Lot_Config', 'Bldg_Type')
  # remove 5 largest houses (before column selection, so it applies to all targets)
  a <- a[a$Gr_Liv_Area <= 4000, , drop = FALSE]
  a <- a[, unique(c(target, z, 'log_price')), drop = FALSE]
  a <- a[
    stats::complete.cases(a[, setdiff(names(a), target), drop = FALSE]),
    ,
    drop = FALSE
  ]
  a <- a[!is.na(a[[target]]), , drop = FALSE]
  if (target == 'Gr_Liv_Area' && nrow(a) != 2925L) {
    warning('Expected 2925 Ames houses for living area; obtained ', nrow(a))
  }
  # Rare frontage factor levels can disappear from observed training cases.
  a$Lot_Config <- as.character(a$Lot_Config)
  a$Lot_Config[a$Lot_Config %in% c('FR2', 'FR3')] <- 'FR2plus'
  a$Lot_Config <- factor(a$Lot_Config)
  a$Bldg_Type <- factor(a$Bldg_Type)
  rownames(a) <- NULL
  attr(a, 'z') <- z
  a
}

ames_formula <- function(response, predictors) {
  stats::reformulate(predictors, response = response)
}

ames_fit <- function(d, target, z, use_y = FALSE) {
  vars <- c(z, if (use_y) 'log_price')
  stats::lm(
    ames_formula(target, vars),
    data = d[!is.na(d[[target]]), , drop = FALSE]
  )
}

ames_fill <- function(d, fit, target, stochastic = FALSE) {
  ix <- which(is.na(d[[target]]))
  if (!length(ix)) {
    return(d)
  }
  mu <- as.numeric(stats::predict(fit, newdata = d[ix, , drop = FALSE]))
  if (stochastic) {
    mu <- mu + stats::rnorm(length(ix), sd = summary(fit)$sigma)
  }
  d[[target]][ix] <- mu
  d
}

ames_mi_train <- function(d, target, z, m, seed) {
  if (!requireNamespace('mice', quietly = TRUE)) {
    stop('Install mice')
  }
  p <- matrix(0L, ncol(d), ncol(d), dimnames = list(names(d), names(d)))
  p[target, c(z, 'log_price')] <- 1L
  method <- setNames(rep('', ncol(d)), names(d))
  method[target] <- 'norm'
  imp <- mice::mice(
    d,
    m = m,
    maxit = 5,
    method = method,
    predictorMatrix = p,
    printFlag = FALSE,
    seed = seed
  )
  lapply(seq_len(m), function(j) mice::complete(imp, j))
}

ames_prediction_draws <- function(
  train_list,
  test_list,
  target,
  price_formula,
  masked_test,
  ytrue,
  probs = c(.05, .95)
) {
  # Mixture of fitted conditional normal price distributions. Includes
  # residual outcome noise and variation across completed datasets/predictions.
  # Does NOT fully integrate regression/imputation model parameter uncertainty.
  n <- nrow(test_list[[1]])
  m <- length(train_list)
  mu <- matrix(NA_real_, n, m)
  sigma <- numeric(m)
  for (j in seq_len(m)) {
    fit <- stats::lm(price_formula, data = train_list[[j]])
    mu[, j] <- as.numeric(stats::predict(fit, newdata = test_list[[j]]))
    sigma[j] <- summary(fit)$sigma
  }
  mean_prediction <- rowMeans(mu)
  # Evaluate both groups of held-out test rows using the same fitted
  # predictive mixture. For observed X, all imputations retain the actual X.
  # For masked X, imputations integrate over uncertainty in X.
  mixq <- function(row, p) {
    means <- mu[row, ]
    spread <- max(sigma)
    lower <- min(means) - 8 * spread
    upper <- max(means) + 8 * spread
    cdf <- function(y) mean(stats::pnorm(y, mean = means, sd = sigma))
    stats::uniroot(function(y) cdf(y) - p, c(lower, upper), tol = 1e-6)$root
  }
  lower <- vapply(seq_len(n), mixq, numeric(1), p = probs[1])
  upper <- vapply(seq_len(n), mixq, numeric(1), p = probs[2])
  summarize_group <- function(ii) {
    errors <- ytrue[ii] - mean_prediction[ii]
    c(
      price_rmse = sqrt(mean(errors^2)),
      price_mae = mean(abs(errors)),
      pi_coverage = mean(ytrue[ii] >= lower[ii] & ytrue[ii] <= upper[ii]),
      pi_width = mean(upper[ii] - lower[ii])
    )
  }
  missing <- summarize_group(which(masked_test))
  observed <- summarize_group(which(!masked_test))
  # Preserve original columns for masked rows; append observed-row metrics.
  c(missing, setNames(observed, paste0(names(observed), '_observed')))
}

# Train/fill predictive FCS for Goals 2 and 3 only.
# Training outcomes are observed except for a small, independently masked
# subset that forces mice(tasks='train') to store a model for log_price.
# All test outcomes are NA during filling. This method predicts Y using
# X and Z; unlike the original Y~X pipelines it also uses Z directly.
ames_fcs_predict <- function(train, test, target, z, masked_test, ytrue,
                             m_train, m_fill, seed, y_train_mask = 0.10,
                             maxit_train = 10L, maxit_fill = 20L) {
  if (!requireNamespace('mice', quietly = TRUE) ||
      !('tasks' %in% names(formals(mice::mice)))) {
    stop('FCS_train_fill requires a mice development version with tasks=train/fill')
  }
  stopifnot(m_train >= 2L, m_fill >= 2L, y_train_mask > 0,
            y_train_mask < 1, all(!is.na(train$log_price)))
  # Factor levels and column order are inherited from the training frame.
  tr <- train
  te <- test
  set.seed(seed)
  mask_y <- sample.int(nrow(tr), max(1L, round(nrow(tr) * y_train_mask)))
  tr$log_price[mask_y] <- NA_real_
  te$log_price <- NA_real_  # strict outcome holdout
  vars <- names(tr)
  pm <- matrix(0L, length(vars), length(vars), dimnames = list(vars, vars))
  pm[target, c(z, 'log_price')] <- 1L
  pm['log_price', c(target, z)] <- 1L
  meth <- setNames(rep('', length(vars)), vars)
  meth[c(target, 'log_price')] <- 'norm'
  trained <- mice::mice(tr, m = m_train, maxit = maxit_train,
                       method = meth, predictorMatrix = pm,
                       tasks = 'train', printFlag = FALSE, seed = seed + 1L)
  if (is.null(trained$models) ||
      !all(c(target, 'log_price') %in% names(trained$models))) {
    stop('FCS train did not store both X and Y models')
  }
  # Workaround: mice.impute.norm does not record model$class (pmm etc. do),
  # so scan.newdata() would refuse to 'fill' with these models.
  for (v in c(target, 'log_price')) {
    trained$models[[v]] <- lapply(trained$models[[v]], function(mod) {
      if (is.null(mod$class)) mod$class <- class(tr[[v]])[1L]
      mod
    })
  }
  # Only X and Y have stored models; complete Z columns are left as is.
  tasks <- setNames(rep('impute', length(vars)), vars)
  tasks[c(target, 'log_price')] <- 'fill'
  filled <- mice::mice(te, m = m_fill, maxit = maxit_fill,
                      method = meth, predictorMatrix = pm,
                      tasks = tasks, models = trained$models,
                      printFlag = FALSE, seed = seed + 2L)
  n <- nrow(te)
  ydraw <- vapply(seq_len(m_fill), function(j) mice::complete(filled, j)$log_price,
                  numeric(n))
  if (anyNA(ydraw)) stop('FCS fill returned missing outcome draws')
  # FCS must never overwrite a genuinely observed test predictor.
  observed <- which(!masked_test)
  if (length(observed) && !all(vapply(seq_len(m_fill), function(j) {
    isTRUE(all.equal(mice::complete(filled, j)[[target]][observed],
                     test[[target]][observed], check.attributes = FALSE))
  }, logical(1)))) stop('FCS fill modified observed test X')
  yhat <- rowMeans(ydraw)
  lo <- apply(ydraw, 1L, stats::quantile, probs = .05, names = FALSE)
  hi <- apply(ydraw, 1L, stats::quantile, probs = .95, names = FALSE)
  score <- function(ii) {
    if (!length(ii)) return(setNames(rep(NA_real_, 4L),
                                      c('price_rmse','price_mae','pi_coverage','pi_width')))
    err <- ytrue[ii] - yhat[ii]
    c(price_rmse = sqrt(mean(err^2)), price_mae = mean(abs(err)),
      pi_coverage = mean(ytrue[ii] >= lo[ii] & ytrue[ii] <= hi[ii]),
      pi_width = mean(hi[ii] - lo[ii]))
  }
  a <- score(which(masked_test))
  b <- score(which(!masked_test))
  c(a, setNames(b, paste0(names(b), '_observed')))
}

ames_slope <- function(d, target) {
  # Descriptive complete-data slope discrepancy, no population bias claim.
  unname(stats::coef(stats::lm(ames_formula('log_price', target), data = d))[2])
}

ames_one <- function(
  target,
  rep,
  m,
  missing_fraction,
  test_fraction,
  seed,
  price_model = 'linear',
  fcs = TRUE,
  m_train = 100L,
  m_fill = 1000L
) {
  d <- ames_preprocess(target)
  z <- attr(d, 'z')
  n <- nrow(d)
  set.seed(seed + rep)
  test_ids <- sample.int(n, round(n * test_fraction))
  truth_train <- d[-test_ids, , drop = FALSE]
  truth_test <- d[test_ids, , drop = FALSE]
  train <- truth_train
  test <- truth_test
  train_ids <- sample.int(nrow(train), round(nrow(train) * missing_fraction))
  test_mask <- sample.int(nrow(test), round(nrow(test) * missing_fraction))
  masked_test <- rep(FALSE, nrow(test))
  masked_test[test_mask] <- TRUE
  train[[target]][train_ids] <- NA_real_
  test[[target]][test_mask] <- NA_real_
  # Prediction model for Goals 2-3: X (linear or quadratic) plus Z.
  # Goal 4 uses ames_slope(), i.e. log_price ~ X only.
  f <- if (price_model == 'linear') {
    ames_formula('log_price', c(target, z))
  } else {
    ames_formula('log_price', c(target, paste0('I(', target, '^2)'), z))
  }
  xfit_z <- ames_fit(train, target, z, FALSE)
  xfit_zy <- ames_fit(train, target, z, TRUE)
  train_mean_z <- ames_fill(train, xfit_z, target)
  train_mean_zy <- ames_fill(train, xfit_zy, target)
  test_mean_z <- ames_fill(test, xfit_z, target)
  set.seed(seed + 10000L * rep + 1L)
  train_draw_z <- ames_fill(train, xfit_z, target, TRUE)
  test_draw_z <- ames_fill(test, xfit_z, target, TRUE)
  train_mi <- ames_mi_train(train, target, z, m, seed + 10000L * rep + 2L)
  set.seed(seed + 10000L * rep + 5L)
  train_draw_zy <- ames_fill(train, xfit_zy, target, TRUE)
  set.seed(seed + 10000L * rep + 3L)
  test_draw_zy <- ames_fill(test, xfit_z, target, TRUE)
  set.seed(seed + 10000L * rep + 4L)
  test_mi <- replicate(
    m,
    ames_fill(test, xfit_z, target, TRUE),
    simplify = FALSE
  )
  pipelines <- list(
    LM_mean_Z = list(train = list(train_mean_z), test = list(test_mean_z)),
    LM_mean_ZY = list(train = list(train_mean_zy), test = list(test_mean_z)),
    LM_draw_Z = list(train = list(train_draw_z), test = list(test_draw_z)),
    LM_draw_ZY = list(train = list(train_draw_zy), test = list(test_draw_zy)),
    LM_MI_ZY = list(train = train_mi, test = test_mi),
    Oracle = list(train = list(truth_train), test = list(truth_test))
  )
  # Reconstruction and descriptive slope diagnostics use training data.
  # Oracle uses true X and the same Y ~ X + Z predictive specification.
  # Reconstruction RMSE is evaluated on artificially hidden training X only.
  full_slope <- ames_slope(truth_train, target)
  out <- lapply(names(pipelines), function(name) {
    p <- pipelines[[name]]
    xhat <- vapply(
      p$train,
      function(t) t[[target]][train_ids],
      numeric(length(train_ids))
    )
    if (is.null(dim(xhat))) {
      xhat <- matrix(xhat, ncol = 1)
    }
    # Per-draw reconstruction loss (not loss of MI mean).
    recon_rmse <- mean(sqrt(colMeans(
      (xhat - truth_train[[target]][train_ids])^2
    )))
    # SD of completed training X, averaged across imputations.
    sd_completed <- mean(vapply(
      p$train,
      function(t) stats::sd(t[[target]]),
      numeric(1)
    ))
    sd_original <- stats::sd(truth_train[[target]])
    variance_ratio <- mean(vapply(
      p$train,
      function(t) {
        stats::var(t[[target]]) / stats::var(truth_train[[target]])
      },
      numeric(1)
    ))
    slope_diff <- mean(vapply(
      p$train,
      function(t) {
        (ames_slope(t, target) - full_slope) / full_slope
      },
      numeric(1)
    ))
    prediction <- ames_prediction_draws(
      p$train,
      p$test,
      target,
      f,
      masked_test,
      truth_test$log_price
    )
    data.frame(
      target = target,
      rep = rep,
      price_model = price_model,
      method = name,
      n_complete = n,
      n_masked_train = length(train_ids),
      n_masked_test = length(test_mask),
      reconstruction_rmse = recon_rmse,
      variance_ratio = variance_ratio,
      sd_original = sd_original,
      sd_completed = sd_completed,
      relative_slope_difference = slope_diff,
      price_rmse = unname(prediction['price_rmse']),
      price_mae = unname(prediction['price_mae']),
      pi_coverage = unname(prediction['pi_coverage']),
      pi_width = unname(prediction['pi_width']),
      price_rmse_observed = unname(prediction['price_rmse_observed']),
      price_mae_observed = unname(prediction['price_mae_observed']),
      pi_coverage_observed = unname(prediction['pi_coverage_observed']),
      pi_width_observed = unname(prediction['pi_width_observed'])
    )
  })
  if (fcs) {
    # Goal 1 reconstruction and Goal 4 slope are not evaluated for this
    # prediction-only procedure. NA prevents misleading comparisons.
    fcs_pred <- ames_fcs_predict(
      train, test, target, z, masked_test, truth_test$log_price,
      m_train = m_train, m_fill = m_fill,
      seed = seed + 10000L * rep + 600L
    )
    out[[length(out) + 1L]] <- data.frame(
      target = target, rep = rep, price_model = price_model,
      method = 'FCS_train_fill', n_complete = n,
      n_masked_train = length(train_ids), n_masked_test = length(test_mask),
      reconstruction_rmse = NA_real_, variance_ratio = NA_real_,
      sd_original = stats::sd(truth_train[[target]]),
      sd_completed = NA_real_, relative_slope_difference = NA_real_,
      price_rmse = unname(fcs_pred['price_rmse']),
      price_mae = unname(fcs_pred['price_mae']),
      pi_coverage = unname(fcs_pred['pi_coverage']),
      pi_width = unname(fcs_pred['pi_width']),
      price_rmse_observed = unname(fcs_pred['price_rmse_observed']),
      price_mae_observed = unname(fcs_pred['price_mae_observed']),
      pi_coverage_observed = unname(fcs_pred['pi_coverage_observed']),
      pi_width_observed = unname(fcs_pred['pi_width_observed'])
    )
  }
  do.call(rbind, out)
}

run_ames_comparison <- function(
  targets = AMES_CANDIDATES,
  nrep = 10L,
  m = 20L,
  missing_fraction = 1 / 3,
  test_fraction = .25,
  seed = 20261022L,
  price_model = 'linear',
  verbose = TRUE,
  fcs = TRUE,
  m_train = 100L,
  m_fill = 1000L
) {
  stopifnot(
    nrep >= 1,
    m >= 2,
    missing_fraction > 0,
    missing_fraction < 1,
    test_fraction > 0,
    test_fraction < 1
  )
  result <- list()
  k <- 0L
  for (target in targets) {
    for (b in seq_len(nrep)) {
      if (verbose && (b == 1L || b %% 10L == 0L)) {
        message(target, ' replication ', b, '/', nrep)
      }
      k <- k + 1L
      result[[k]] <- ames_one(
        target,
        b,
        m,
        missing_fraction,
        test_fraction,
        seed,
        price_model,
        fcs = fcs, m_train = m_train, m_fill = m_fill
      )
    }
  }
  ans <- do.call(rbind, result)
  rownames(ans) <- NULL
  ans
}

summary_table <- function(results) {
  metrics <- c(
    'reconstruction_rmse',
    'variance_ratio',
    'sd_original',
    'sd_completed',
    'relative_slope_difference',
    'price_rmse',
    'price_mae',
    'pi_coverage',
    'pi_width',
    'price_rmse_observed',
    'price_mae_observed',
    'pi_coverage_observed',
    'pi_width_observed'
  )
  stats::aggregate(
    results[, metrics],
    by = list(target = results$target, method = results$method),
    FUN = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  )
}

# Preliminary measure of predictive strength on the fixed complete-case data.
# R2 of Y~X and incremental R2 over Y~Z (descriptive, not cross-validated).
predictor_strength <- function(targets = AMES_CANDIDATES) {
  do.call(
    rbind,
    lapply(targets, function(target) {
      d <- ames_preprocess(target)
      z <- attr(d, 'z')
      r2x <- summary(stats::lm(
        ames_formula('log_price', target),
        data = d
      ))$r.squared
      r2z <- summary(stats::lm(
        ames_formula('log_price', z),
        data = d
      ))$r.squared
      r2zx <- summary(stats::lm(
        ames_formula('log_price', c(z, target)),
        data = d
      ))$r.squared
      data.frame(
        target = target,
        n = nrow(d),
        r2_x = r2x,
        r2_z = r2z,
        r2_z_plus_x = r2zx,
        incremental_r2 = r2zx - r2z
      )
    })
  )
}

# Convenient display of the two coverage groups side by side.
coverage_by_status <- function(results) {
  s <- summary_table(results)
  s[, c(
    'target',
    'method',
    'pi_coverage_observed',
    'pi_coverage',
    'pi_width_observed',
    'pi_width'
  )]
}
