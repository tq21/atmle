r_loss_cate <- function(S,
                        A,
                        Y,
                        theta,
                        Q_bar,
                        g1W,
                        Pi1WA,
                        phi,
                        beta) {
  return(Y$sub(theta)$sub(A$sub(g1W)$mul(phi$matmul(beta)))$pow(2)$mean())
}
