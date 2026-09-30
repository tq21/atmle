# Fixed initial HAL working spaces. Always retain the intercept and any
# explicitly unpenalized treatment column, even if its fitted coefficient is 0.
.atmle_working_models <- function(fit, n_lambda, source = FALSE) {
  is_glm <- identical(fit$fit_type, "glm")
  if (is_glm) {
    lambda <- NA_real_
    beta_cv <- as.numeric(fit$fit$coefficients)
    beta_cv[is.na(beta_cv)] <- 0
  } else {
    lambda_cv <- fit$fit$lambda.min
    lambda <- fit$fit$lambda[fit$fit$lambda <= lambda_cv]
    lambda <- utils::head(lambda, n_lambda)
    beta_cv <- as.numeric(stats::coef(fit$fit, s = lambda_cv))
  }
  forced <- fit$forced_columns
  selected <- function(b) {
    if (is_glm) seq_along(b) else sort(unique(c(1L, which(b != 0), forced)))
  }
  cv_columns <- selected(beta_cv)
  design <- cbind(1, fit$phi)
  if (source) {
    design1 <- cbind(1, fit$phi_W1)
    design0 <- cbind(1, fit$phi_W0)
  }
  lapply(seq_along(lambda), function(j) {
    columns <- if (is_glm) seq_along(beta_cv) else {
      selected(as.numeric(stats::coef(fit$fit, s = lambda[j])))
    }
    model <- list(idx = j, lambda = lambda[j], beta = beta_cv[columns],
                  pseudo_outcome = fit$pseudo_outcome,
                  pseudo_weights = fit$pseudo_weights)
    if (source) {
      model$phi_WA <- design[, columns, drop = FALSE]
      model$phi_W1 <- design1[, columns, drop = FALSE]
      model$phi_W0 <- design0[, columns, drop = FALSE]
      model$cate_WA <- as.numeric(design[, cv_columns, drop = FALSE] %*% beta_cv[cv_columns])
      model$cate_W1 <- as.numeric(design1[, cv_columns, drop = FALSE] %*% beta_cv[cv_columns])
      model$cate_W0 <- as.numeric(design0[, cv_columns, drop = FALSE] %*% beta_cv[cv_columns])
    } else model$phi_W <- design[, columns, drop = FALSE]
    model
  })
}

.atmle_validate_inference <- function(alpha, weight_cap) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) stop("alpha must lie strictly between 0 and 1.")
  if (!is.numeric(weight_cap) || length(weight_cap) != 1L ||
      !is.finite(weight_cap) || weight_cap < 0 || weight_cap > 1) {
    stop("weight_cap must lie in [0, 1].")
  }
}

# Target both source-bias candidates using the same fitted nuisances and data.
.atmle_target_pair <- function(fit, source_fits, verbose = FALSE) {
  workers <- list()
  components <- list()
  shared_fields <- c("controls_only", "beta_target_method", "target_gwt", "settings",
                     "A", "Y", "S", "Delta", "weights", "g_bar", "g_bar0", "theta",
                     "Pi_bar", "Pi", "Q_bar", "tau_A")
  for (candidate in c("forced", "unforced")) {
    worker <- atmle_ate_fusion$new(fit$data, fit$S_node, fit$W_nodes,
                                   fit$A_node, fit$Y_node, fit$family, fit$n_folds)
    for (field in shared_fields) worker[[field]] <- fit[[field]]
    worker$settings$primary <- candidate
    worker$tau_S <- source_fits[[candidate]]
    if (candidate == "forced") {
      worker$tau_A_star <- worker$target_beta_A_seq(worker$settings$n_lambda, FALSE)
      worker$tau_A_star_avg_over_S1 <- worker$target_beta_A_seq(worker$settings$n_lambda, TRUE)
    } else {
      worker$tau_A_star <- workers$forced$tau_A_star
      worker$tau_A_star_avg_over_S1 <- workers$forced$tau_A_star_avg_over_S1
    }
    worker$tau_S_star <- worker$target_Pi_beta_S(worker$settings$n_lambda, FALSE,
                                                worker$settings$max_iter, verbose)
    worker$tau_S_star_avg_over_S1 <- worker$target_Pi_beta_S(worker$settings$n_lambda, TRUE,
                                                            worker$settings$max_iter, verbose)
    worker$Pi_star <- lapply(worker$tau_S_star, `[[`, "Pi_star")
    worker$Pi_star_avg_over_S1 <- lapply(worker$tau_S_star_avg_over_S1, `[[`, "Pi_star")
    worker$inference(worker$settings$alpha)
    workers[[candidate]] <- worker
    components[[candidate]] <- list(results = worker$results, influence = worker$influence)
  }
  list(workers = workers, components = components)
}

# All variances below are variances of estimators (IF variance divided by n).
.atmle_floor_pair <- function(unforced, forced, du, df, weight_cap = 1) {
  if (length(du) != length(df) || length(du) < 2L ||
      any(!is.finite(c(unforced, forced, du, df)))) {
    stop("The combination requires finite estimates and aligned influence curves.")
  }
  n <- length(du)
  vu <- stats::var(du) / n
  vf <- stats::var(df) / n
  covariance <- stats::cov(du, df) / n
  vd <- stats::var(du - df) / n
  gap <- unforced - forced
  bias_squared <- max(gap^2, vd)
  tiny <- 1e-10 * max(vu, vf, .Machine$double.xmin)
  weight <- if (vd <= tiny) 0 else {
    min(weight_cap, max(0, (vf - covariance) / (vd + bias_squared)))
  }
  influence <- (1 - weight) * df + weight * du
  list(psi = forced + weight * gap, influence = influence,
       se = sqrt(stats::var(influence) / n), weight = weight,
       forced_variance = vf, unforced_variance = vu,
       component_covariance = covariance, difference_variance = vd,
       bias_squared = bias_squared)
}

.atmle_influence <- function(fit) {
  get <- function(a, s) {
    unlist(lapply(a, function(aa) lapply(s, function(ss) aa$eic - ss$eic)),
           recursive = FALSE)
  }
  c(get(fit$tau_A_star, fit$tau_S_star),
    get(fit$tau_A_star_avg_over_S1, fit$tau_S_star_avg_over_S1))
}

.atmle_finish_inference <- function(fit, alpha) {
  .atmle_validate_inference(alpha, fit$settings$weight_cap)
  fit$influence <- .atmle_influence(fit)
  primary <- fit$settings$primary
  fit$results$estimator <- paste0("A-TMLE (", primary, " A)")
  fits <- c(lapply(fit$tau_A_star, function(x) fit$tau_S_star),
            lapply(fit$tau_A_star_avg_over_S1, function(x) fit$tau_S_star_avg_over_S1))
  source_ok <- unlist(lapply(fits, function(x) vapply(x, function(y) isTRUE(y$converged), logical(1))))
  finite <- vapply(fit$influence, function(x) all(is.finite(x)), logical(1))
  fit$results$converged <- source_ok & finite & is.finite(fit$results$psi)
  for (field in c("unforced_weight", "forced_variance", "unforced_variance",
                  "component_covariance", "difference_variance", "bias_squared")) {
    fit$results[[field]] <- NA_real_
  }
  if (!is.null(fit$components)) {
    report <- .atmle_combine_report(fit$components, primary, alpha, fit$settings$weight_cap)
    fit$results <- report$results
    fit$influence <- report$influence
  }
  invisible(fit)
}

.atmle_combine_report <- function(components, primary, alpha, weight_cap) {
  results <- components[[primary]]$results
  influence <- components[[primary]]$influence
  z <- stats::qnorm(1 - alpha/2)
  results$alpha <- alpha
  results$lower <- results$psi - z * results$se
  results$upper <- results$psi + z * results$se
  u <- components$unforced
  f <- components$forced
  for (population in unique(results$param)) {
    # One additional row per population, using each candidate's CV-selected
    # working models. Extra undersmoothing rows are retained separately.
    iu <- which(u$results$param == population & u$results$tau_A_idx == 1L &
                  u$results$tau_S_idx == 1L)
    iff <- which(f$results$param == population & f$results$tau_A_idx == 1L &
                   f$results$tau_S_idx == 1L)
    if (length(iu) != 1L || length(iff) != 1L) stop("Missing CV-selected candidate pair.")
    combo <- .atmle_floor_pair(u$results$psi[iu], f$results$psi[iff],
                               u$influence[[iu]], f$influence[[iff]], weight_cap)
    row <- f$results[iff, , drop = FALSE]
    row$estimator <- "Variance-floor convex combination"
    row$tau_A_lambda <- row$tau_S_lambda <- NA_real_
    row$psi_tilde <- (1-combo$weight)*f$results$psi_tilde[iff] + combo$weight*u$results$psi_tilde[iu]
    row$psi_pound <- (1-combo$weight)*f$results$psi_pound[iff] + combo$weight*u$results$psi_pound[iu]
    row$psi <- combo$psi
    row$se <- combo$se
    row$alpha <- alpha
    row$lower <- combo$psi - z*combo$se
    row$upper <- combo$psi + z*combo$se
    row$unforced_weight <- combo$weight
    row$converged <- u$results$converged[iu] && f$results$converged[iff]
    for (field in c("forced_variance", "unforced_variance", "component_covariance",
                    "difference_variance", "bias_squared")) row[[field]] <- combo[[field]]
    results <- rbind(results, row)
    influence[[length(influence) + 1L]] <- combo$influence
  }
  rownames(results) <- NULL
  list(results = results, influence = influence)
}
