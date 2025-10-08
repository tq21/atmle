#' @title Nonparametric TMLE for Psi (std. over pooled-covariates) and
#' Psi_2 (std. over RCT-covariates)
#'
#' @importFrom purrr map
#'
#' @export
nonparametric <- function(data,
                          S,
                          W,
                          A,
                          Y,
                          family,
                          Pi_method = "glm",
                          g_method = "glm",
                          Q_method = "glm",
                          Q_pooling = TRUE,
                          v_folds = 5,
                          verbose = FALSE,
                          browse = FALSE) {

  if (browse) browser()

  # define nodes ---------------------------------------------------------------
  S <- data[[S]]
  W <- data[, W, drop = FALSE]
  A <- data[[A]]
  Y <- data[[Y]]
  n <- nrow(data) # sample size
  g_bounds <- c(5/sqrt(n)/log(n), 1-5/sqrt(n)/log(n))
  Pi_bounds <- c(5/sqrt(n)/log(n), 1)

  # cross fitting schemes
  if (family == "gaussian") {
    cv_strata <- paste0(S, "-", A)
  } else if (family == "binomial") {
    cv_strata <- paste0(S, "-", A, "-", Y)
  }
  suppressWarnings({
    folds <- make_folds(n = n, V = v_folds,
                        strata_ids = as.integer(factor(cv_strata)))
  })
  foldid <- folds2foldvec(folds)
  foldid_S1 <- foldid[S == 1]
  folds_S1 <- map(seq(v_folds), function(v) {
    fold_from_foldvec(v = v, folds = foldid_S1)
  })

  # estimate nuisance parts
  if (verbose) print("learning E(Y|S,W,A)")
  Q <- learn_QSWA(S = S,
                  W = W,
                  A = A,
                  Y = Y,
                  folds = folds,
                  folds_S1 = folds_S1,
                  family = family,
                  method = Q_method,
                  pooling = Q_pooling)

  if (verbose) print("learning P(A=1|S,W)")
  g11W <- learn_g11W(S = S,
                     W = W,
                     A = A,
                     method = g_method,
                     g_bounds = g_bounds)

  if (verbose) print("learning P(S=1|W)")
  Pi <- learn_SW(S = S,
                 W = W,
                 folds = folds,
                 method = Pi_method,
                 Pi_bounds = Pi_bounds)

  # target Q (pooled W)
  Q_star <- target_Q(S = S,
                     W = W,
                     A = A,
                     Y = Y,
                     Pi = Pi,
                     g11W = g11W,
                     Q = Q)

  # target Q (RCT W)
  pS <- mean(S)
  Q_star_rct_W <- target_Q_rct_W(S = S,
                                 W = W,
                                 A = A,
                                 Y = Y,
                                 pS = pS,
                                 g11W = g11W,
                                 Q = Q)

  # parameter that avg. over pooled W
  psi_pooled_W <- mean(Q_star$Q1W1-Q_star$Q1W0)
  eic_pooled_W <- get_np_eic_pooled_W(Q = Q_star,
                                      Pi = Pi,
                                      g11W = g11W,
                                      S = S,
                                      A = A,
                                      Y = Y,
                                      psi = psi_pooled_W)
  se_pooled_W <- sqrt(var(eic_pooled_W, na.rm = TRUE)/n)
  lower_pooled_W <- psi_pooled_W+qnorm(0.025)*se_pooled_W
  upper_pooled_W <- psi_pooled_W+qnorm(0.975)*se_pooled_W

  # parameter that avg. over RCT W
  psi_rct_W <- weighted.mean(Q_star_rct_W$Q1W1[S==1]-Q_star_rct_W$Q1W0[S==1],
                             w = (S/pS)[S==1])
  eic_rct_W <- get_np_eic_rct_W(Q = Q_star_rct_W,
                                pS = pS,
                                g11W = g11W,
                                S = S,
                                A = A,
                                Y = Y,
                                psi = psi_rct_W)
  se_rct_W <- sqrt(var(eic_rct_W, na.rm = TRUE)/n)
  lower_rct_W <- psi_rct_W+qnorm(0.025)*se_rct_W
  upper_rct_W <- psi_rct_W+qnorm(0.975)*se_rct_W

  return(list(psi_pooled_W = psi_pooled_W,
              lower_pooled_W = lower_pooled_W,
              upper_pooled_W = upper_pooled_W,
              eic_pooled_W = eic_pooled_W,
              psi_rct_W = psi_rct_W,
              lower_rct_W = lower_rct_W,
              upper_rct_W = upper_rct_W,
              eic_rct_W = eic_rct_W))
}
