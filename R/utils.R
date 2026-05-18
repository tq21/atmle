.bound <- function(x, bounds) {
  pmin(pmax(x, bounds[1]), bounds[2])
}

mat_inverse <- function(mat,
                        small_diag = 1e-3,
                        fall_back_method = "svd_pseudo_inv") {
  tryCatch({
    solve(mat)
  }, error = function(e) {
    if (fall_back_method == "svd_pseudo_inv") {
      svd_pseudo_inv(mat)
    } else if (fall_back_method == "ridge") {
      solve(mat + diag(small_diag, nrow(mat), ncol(mat)))
    } else {
      stop("Unknown fall_back_method specified.")
    }
  })
}

svd_pseudo_inv <- function(mat, tol = 1e-3) {
  svd_res <- svd(mat)
  D_inv <- ifelse(svd_res$d > tol, 1 / svd_res$d, 0)
  return(svd_res$v %*% diag(D_inv) %*% t(svd_res$u))
}
