r_loss_care_with_Pi <- function(S,
                                A,
                                Y,
                                theta,
                                Q_bar,
                                g1W,
                                Pi1WA,
                                phi,
                                beta,
                                epsilon,
                                ones,
                                eps) {

  r_loss <- Y$sub(Q_bar)$sub(S$sub(Pi1WA)$mul(phi$matmul(beta)))$pow(2)
  H <- A$div(g1W)$sub(ones$sub(A)$div(ones$sub(g1W)))$mul(phi$matmul(beta))
  targeted_loss_Pi <- S*torch_log(torch_sigmoid(Pi1WA+epsilon*H)+eps)+
    (1-S)*torch_log(1-torch_sigmoid(Pi1WA+epsilon*H)+eps)

  return(r_loss$sub(targeted_loss_Pi)$mean())
}
