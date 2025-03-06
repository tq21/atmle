grad_proj <- function(W,
                      A,
                      Y,
                      theta,
                      g1W,
                      phi_W,
                      beta) {

  # obtain score matrix
  cate_pred <- as.numeric(phi_W%*%beta)
  score_mat <- (Y-theta-(A-g1W)*cate_pred)*(A-g1W)*phi_W

  # compute canonical gradient
  eic_np <- eic_ate(QW1 = theta+(1-g1W)*cate_pred,
                    QW0 = theta-g1W*cate_pred,
                    psi = mean(cate_pred),
                    A = A,
                    g1W = g1W,
                    Y = Y,
                    QWA = theta+(A-g1W)*cate_pred)

  # project gradient onto the score space (regularized)
  proj_fit <- glmnet(x = score_mat,
                     y = eic_np,
                     intercept = FALSE,
                     lambda = 1e-5,
                     alpha = 1)
  proj_eic <- as.numeric(predict(proj_fit, newx = score_mat, s = 1e-5))

  return(proj_eic)
}
