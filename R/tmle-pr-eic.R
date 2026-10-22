# Internal influence curves retained from the original TMLE-PR implementation.
get_np_eic_pooled_W <- function(Q,
                                Pi,
                                g11W,
                                S,
                                A,
                                Y,
                                psi,
                                g_bound = NULL) {
  if (is.null(g_bound)) {
    g_bound <- 5/sqrt(length(Y))/log(length(Y))
  }
  g01W <- 1-g11W
  denom_A1 <- .bound(Pi*g11W, c(g_bound, 1))
  denom_A0 <- .bound(Pi*g01W, c(g_bound, 1))
  W_comp <- Q$Q1W1-Q$Q1W0-psi
  Q_comp <- S*(A/denom_A1-(1-A)/denom_A0)*(Y-Q$Q1WA)
  return(W_comp+Q_comp)
}

get_np_eic_rct_W <- function(Q,
                             pS,
                             g11W,
                             S,
                             A,
                             Y,
                             psi,
                             g_bound = NULL) {
  if (is.null(g_bound)) {
    g_bound <- 5/sqrt(length(Y))/log(length(Y))
  }
  g01W <- 1-g11W
  denom_A1 <- .bound(pS*g11W, c(g_bound, 1))
  denom_A0 <- .bound(pS*g01W, c(g_bound, 1))
  W_comp <- (S/pS)*(Q$Q1W1-Q$Q1W0-psi)
  Q_comp <- S*(A/denom_A1-(1-A)/denom_A0)*(Y-Q$Q1WA)
  return(W_comp+Q_comp)
}
