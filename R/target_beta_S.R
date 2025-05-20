target_beta_S <- function(S,
                          W,
                          A,
                          Y,
                          g1W,
                          Pi,
                          theta_WA,
                          tau_S,
                          weights,
                          controls_only,
                          target_method,
                          browse = FALSE) {
  if (browse) browser()

  # targeting
  if (target_method == "relaxed") {
    # relaxed-fit targeting
    relaxed_fit <- glm(tau_S$pseudo_outcome ~ -1+.,
                       family = "gaussian",
                       data = as.data.frame(tau_S$phi_WA),
                       weights = tau_S$pseudo_weights)
    tau_S$beta <- as.numeric(coef(relaxed_fit))
    na_idx <- which(is.na(tau_S$beta))
    if (length(na_idx) > 0) {
      tau_S$beta <- tau_S$beta[!is.na(tau_S$beta)]
      tau_S$phi_WA <- tau_S$phi_WA[, -na_idx, drop = FALSE]
      tau_S$phi_W0 <- tau_S$phi_W0[, -na_idx, drop = FALSE]
      tau_S$phi_W1 <- tau_S$phi_W1[, -na_idx, drop = FALSE]
    }
  } else if (target_method == "oneshot") {
    # TODO: implement target in a sequence of WMs
    IM <- t(phi_WA) %*% diag(Pi$pred*(1-Pi$pred)) %*% phi_WA / length(Y)
    IM_inv <- tryCatch({
      solve(IM)
    }, error = function(e) {
      # TODO: add helpful message if verbose
      if (eic_method == "svd_pseudo_inv") {
        svd_pseudo_inv(IM)
      } else if (eic_method == "diag") {
        solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
      } else {
        stop("Unknown eic_method specified.")
      }
    })
    beta <- as.vector(tau_A_obj$beta)
    clever_cov <- as.vector(IM_inv %*% colMeans(phi_W))
    H <- (A-g1W)*as.vector(phi_W %*% clever_cov)
    tau <- as.vector(phi_W %*% beta)
    R <- Y-theta_WA-(A-g1W)*tau
    epsilon <- sum(H*R)/sum(H*H)
    beta <- beta+epsilon*clever_cov
  }

  tau_S$cate_WA <- as.vector(tau_S$phi_WA %*% tau_S$beta)
  tau_S$cate_W0 <- as.vector(tau_S$phi_W0 %*% tau_S$beta)
  tau_S$cate_W1 <- as.vector(tau_S$phi_W1 %*% tau_S$beta)

  return(tau_S)
}
