target_beta_A <- function(A,
                          Y,
                          g1W,
                          theta_W,
                          tau_A,
                          target_method,
                          eic_method,
                          browse = FALSE) {
  if (browse) browser()

  # targeting
  if (target_method == "relaxed") {
    # relaxed-fit targeting
    relaxed_fit <- glm(tau_A$pseudo_outcome ~ -1+.,
                       family = "gaussian",
                       data = data.frame(tau_A$phi_train),
                       weights = tau_A$pseudo_weights)
    tau_A$beta <- as.numeric(coef(relaxed_fit))
    na_idx <- which(is.na(tau_A$beta))
    if (length(na_idx) > 0) {
      tau_A$beta <- tau_A$beta[!is.na(tau_A$beta)]
      tau_A$phi_train <- tau_A$phi_train[, -na_idx, drop = FALSE]
      tau_A$phi_W <- tau_A$phi_W[, -na_idx, drop = FALSE]
    }
  } else if (target_method == "oneshot") {
    # TODO: implement target in a sequence of WMs
    IM <- t(tau_A$phi_W) %*% diag(g1W*(1-g1W)) %*% tau_A$phi_W / length(Y)
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
    clever_cov <- as.vector(IM_inv %*% colMeans(tau_A$phi_W))
    H <- (A-g1W)*as.vector(tau_A$phi_W %*% clever_cov)
    R <- Y-theta_W-(A-g1W)*tau_A$cate_W
    epsilon <- sum(H*R)/sum(H*H)
    tau_A$beta <- tau_A$beta+epsilon*clever_cov
  }

  tau_A$cate_W <- as.vector(tau_A$phi_W %*% tau_A$beta)

  return(tau_A)
}
