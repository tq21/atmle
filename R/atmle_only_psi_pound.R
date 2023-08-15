# max_iter=0
# eps=1e-5
# tmle=TRUE
# psi_tilde_adaptive=TRUE
#
# `%+%` <- function(a, b) paste0(a, b)
#
# atmle <- function(max_iter=10,
#                   eps=1e-5,
#                   tmle=FALSE,
#                   tau_method = "parametric",
#                   psi_pound_method="adaptive",
#                   psi_tilde_method="tmle",
#                   g_rct=0.67) {
#   S <- data$S
#   W <- as.matrix(data[, c("W1", "W2", "W3", "W4")])
#   A <- data$A
#   Y <- data$Y
#
#   n <- nrow(data)
#
#   # estimate bias psi_pound ----------------------------------------------------
#   # learn relevant parts
#   print("learning E(Y|W)")
#   theta_pred <- learn_theta(W, A, Y)
#   print("learning P(S=1|W,A)")
#   Pi_pred <- learn_Pi(S, W, A)
#   print("learning P(A|W)")
#   g_pred <- learn_g(S, W, A, g_rct)
#   print("learning E(Y|S,W,A)")
#   tau_pred <- learn_tau(S, W, A, Y, Pi_pred$pred, theta_pred)
#
#   psi_pound_est <- NULL
#   psi_pound_eic <- NULL
#
#   print("learning working model: tau(W,A)=E(Y|S=1,W,A)-E(Y|S=0,W,A)")
#   if (tau_method == "parametric" | tau_method == "semiparametric") {
#     tau_pred <- learn_tau_parametric(S, W, A, Y)
#   } else if (tau_method = "adaptive") {
#     tau_pred <- learn_tau(S, W, A, Y, Pi_pred$pred, theta_pred)
#   }
#
#
#
#   psi_pound_est <- mean((1-Pi_pred$A0)*tau_pred$A0-(1-Pi_pred$A1)*tau_pred$A1)
#   psi_pound_eic <- get_eic_psi_pound(Pi_pred, tau_pred, psi_pound_est, g_pred, theta_pred, S, A, Y, n)
#
#
#   if (psi_pound_method == "parametric") {
#     # use relaxed HAL to estimate psi_pound
#     tau_pred <- learn_tau_parametric(S, W, A, Y)
#
#
#   } else if (psi_pound_method == "adaptive") {
#     # use adaptive + relaxed HAL plug-in to estimate psi_pound
#     psi_pound_est <- mean((1-Pi_pred$A0)*tau_pred$A0-(1-Pi_pred$A1)*tau_pred$A1)
#     psi_pound_eic <- get_eic_psi_pound(Pi_pred, tau_pred, psi_pound_est, g_pred, theta_pred, S, A, Y, n)
#   } else if (psi_pound_method == "a-tmle") {
#     # use adaptive + TMLE to estimate psi_pound
#     # TODO: perform TMLE update of Pi
#     cur_iter <- 0
#     tol <- Inf
#
#     while (cur_iter <= max_iter && tol > eps) {
#       print(cur_iter)
#       # update
#       Pi_star <- Pi_tmle(S, W, A, tau_pred, Pi_pred)
#
#       # re-learn tau using Pi_star
#       tau_star <- learn_tau(S, W, A, Y, Pi_star$pred, theta_pred)
#
#       # current EIC
#       cur_eic <- get_eic_Pi(g_pred, tau_star, Pi_star, S, A)
#       tol <- abs(mean(cur_eic) - 0)
#
#       cur_iter <- cur_iter + 1
#       print(tol)
#   }
#
#     psi_pound_est <- mean((1-Pi_star$A0)*tau_star$A0-(1-Pi_star$A1)*tau_star$A1)
#     psi_pound_eic <- get_eic_psi_pound(Pi_star, tau_star, g_pred, theta_pred, psi_pound_est, S, A, Y, n)
#
#   } else {
#     # no targeting, just plug-in
#     psi_pound_est <- mean((1-Pi_pred$A0)*tau_pred$A0-(1-Pi_pred$A1)*tau_pred$A1)
#     psi_pound_eic <- get_eic_psi_pound(Pi_pred, tau_pred, psi_pound_est, g_pred, theta_pred, S, A, Y, n)
#   }
#
#   # estimate pooled ATE psi_tilde ----------------------------------------------
#   psi_tilde_est <- NULL
#   psi_tilde_eic <- NULL
#
#   if (psi_tilde_method == "adaptive") {
#     # use A-TMLE to estimate psi_tilde
#     theta_tilde_pred <- learn_theta_tilde(W, Y)
#     psi_tilde <- learn_psi_tilde(W, A, Y, g_pred, theta_tilde_pred)
#
#     psi_tilde_est <- mean(psi_tilde$pred)
#     psi_tilde_eic <- get_eic_psi_tilde(psi_tilde, g_pred, theta_pred, Y, A)
#
#   } else if (psi_tilde_method == "HAL-MLE") {
#     # use relaxed HAL as a plug in to estimate psi_tilde
#     Q_pred <- learn_Q(W, A, Y)
#
#     # TODO: sandwich estimator
#
#   } else if (psi_tilde_method == "tmle") {
#     # use TMLE to estimate psi_tilde
#     Q_pred <- learn_Q(W, A, Y)
#     Q_star <- tmle(Y = Y, A = A, W = W, g1W = g_pred,
#                    Q = as.matrix(data.frame(Q_pred$A1, Q_pred$A0)),
#                    family = "gaussian")
#     psi_tilde_est <- Q_star$estimates$ATE$psi
#     psi_tilde_eic <- Q_star$estimates$IC$IC.ATE
#   }
#
#   # estimate psi ---------------------------------------------------------------
#   psi_est <- psi_tilde_est - psi_pound_est
#   psi_eic <- psi_tilde_eic - psi_pound_eic
#   psi_pound_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
#   psi_ci <- as.character(psi_est-1.96*psi_pound_se) %+% ", " %+% as.character(psi_est+1.96*psi_pound_se)
#
#   print("point estimate: " %+% psi_est %+% ", 95% CI: " %+% psi_ci)
# }
#
# atmle(max_iter=0,eps=1e-5,tmle=FALSE,psi_pound_method="adaptive",psi_tilde_method="tmle")
# atmle(max_iter=0,eps=1e-5,tmle=FALSE,psi_pound_method="adaptive",psi_tilde_method="adaptive")
#
#
# atmle(max_iter=0,eps=1e-5,external_trx=TRUE,tmle=TRUE,psi_tilde_adaptive=FALSE)
# res <- atmle(max_iter=0,eps=1e-5,external_trx=TRUE,tmle=TRUE,psi_tilde_adaptive=TRUE)


