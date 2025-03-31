get_beta_h <- function(x_basis,
                       g1W,
                       eic_method = "svd_pseudo_inv") {
  n <- nrow(x_basis)
  IM <- t(x_basis) %*% diag(g1W*(1-g1W)) %*% x_basis / n
  if (dim(x_basis)[2] == 1) {
    IM_inv <- solve(IM)
  } else {
    if (eic_method == "svd_pseudo_inv") {
      # SVD-based pseudo-inverse
      IM_inv <- svd_pseudo_inv(IM)
    } else if (eic_method == "diag") {
      IM_inv <- solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
    }
  }
  beta_h <- as.vector(IM_inv %*% colMeans(x_basis))
  return(beta_h)
}
