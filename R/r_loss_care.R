r_loss_care <- function(S,
                        A,
                        Y,
                        theta,
                        Q_bar,
                        g1W,
                        Pi1WA,
                        phi_W,
                        phi_WA,
                        beta,
                        alpha) {
  return(Y$sub(Q_bar)$sub(S$sub(Pi1WA)$mul(phi_WA$matmul(alpha)))$pow(2)$mean())
}
