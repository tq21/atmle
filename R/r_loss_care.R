r_loss_care <- function(S,
                        A,
                        Y,
                        theta,
                        Q_bar,
                        g1W,
                        Pi1WA,
                        phi,
                        beta) {
  return(Y$sub(Q_bar)$sub(S$sub(Pi1WA)$mul(phi$matmul(beta)))$pow(2)$mean())
}
