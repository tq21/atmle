library(hal9001)
library(data.table)
library(tmle)
library(glmnet)

# function to learn theta(W,A)=E(Y|W,A)
learn_theta <- function(W, A, Y, method = "lasso") {
  pred <- NULL
  if (method == "lasso") {
    foldid = sample(rep(seq(10), length = length(A)))
    fit <- cv.glmnet(x = as.matrix(data.frame(W, A = A)), y = Y, keep = TRUE, alpha = 1, foldid = foldid)
    y_lambda_min <- fit$lambda.min
    #pred <- as.numeric(predict(fit, newdata = as.matrix(data.frame(W, A = A))))
    pred <- fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == y_lambda_min]
  } else if (method == "HAL") {
    fit <- fit_hal(X = as.matrix(data.frame(W, A = A)), Y = Y, family = "gaussian", smoothness_orders = 0)
    pred <- as.numeric(predict(fit, new_data = as.matrix(data.frame(W, A = A))))
  } else if (method == "glm") {
    fit <- lm(Y ~ ., data = data.frame(W, A = A))
    pred <- as.numeric(predict(fit))
  }

  return(pred)
}

# function to learn theta_tilde(W)=E(Y|W)
learn_theta_tilde <- function(W, Y, method = "lasso") {
  pred <- NULL
  if (method == "lasso") {
    foldid = sample(rep(seq(5), length = length(Y)))
    fit <- cv.glmnet(x = as.matrix(data.frame(W)), y = Y, foldid = foldid, keep = TRUE, alpha = 1)
    y_lambda_min <- fit$lambda.min
    pred <- fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == y_lambda_min]
  } else if (method == "HAL") {
    fit <- fit_hal(X = as.matrix(data.frame(W)), Y = Y, family = "gaussian", smoothness_orders = 0)
    pred <- as.numeric(predict(fit, new_data = as.matrix(data.frame(W))))
  } else if (method == "glm") {
    fit <- lm(Y ~ ., data = data.frame(W))
    pred <- as.numeric(predict(fit))
  }

  return(pred)
}

# function to learn Pi(1|W,A)=P(S=1|W,A)
learn_Pi <- function(S, W, A, method = "lasso") {
  pred <- NULL
  A1 <- NULL
  A0 <- NULL
  if (method == "lasso") {
    foldid = sample(rep(seq(5), length = length(A)))
    fit <- cv.glmnet(x = as.matrix(data.table(W, A = A)), y = S, foldid = foldid, keep = TRUE, alpha = 1, family = "binomial")
    y_lambda_min <- fit$lambda.min
    #pred <- plogis(fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == y_lambda_min])
    pred <- as.numeric(predict(fit, newx = as.matrix(data.table(W, A = A)), s = y_lambda_min, type = "response"))
    A1 <- as.numeric(predict(fit, newx = as.matrix(data.table(W, A = 1)), s = y_lambda_min, type = "response"))
    A0 <- as.numeric(predict(fit, newx = as.matrix(data.table(W, A = 0)), s = y_lambda_min, type = "response"))
  } else if (method == "HAL") {
    fit <- fit_hal(X = as.matrix(data.frame(W, A = A)), Y = S, family = "binomial", smoothness_orders = 0)
    pred <- as.numeric(predict(fit, new_data = as.matrix(data.frame(W, A = A))))
    A1 <- as.numeric(predict(fit, new_data = as.matrix(data.frame(W, A = 1))), type = "response")
    A0 <- as.numeric(predict(fit, new_data = as.matrix(data.frame(W, A = 0))), type = "response")
  } else if (method == "glm") {
    fit <- glm(S ~ ., data = data.frame(W, A = A), family = "binomial")
    pred <- as.numeric(predict(fit, type = "response"))
    A1 <- as.numeric(predict(fit, newdata = data.frame(W, A = 1), type = "response"))
    A0 <- as.numeric(predict(fit, newdata = data.frame(W, A = 0), type = "response"))
  }

  return(list(pred = pred,
              A1 = A1,
              A0 = A0))
}

# function to learn tau, R-learner
learn_tau <- function(S, W, A, Y, Pi, theta, method = "lasso") {
  weights <- 1
  pseudo_outcome_A0 <- (Y[A == 0]-theta[A == 0])/(S[A == 0]-Pi$pred[A == 0])
  pseudo_weights_A0 <- (S[A == 0]-Pi$pred[A == 0])^2*weights
  pseudo_outcome_A1 <- (Y[A == 1]-theta[A == 1])/(S[A == 1]-Pi$pred[A == 1])
  pseudo_weights_A1 <- (S[A == 1]-Pi$pred[A == 1])^2*weights

  pred <- NULL
  A1 <- NULL
  A0 <- NULL
  x_basis <- NULL
  x_basis_A1 <- NULL
  x_basis_A0 <- NULL
  if (method == "lasso") {
    # A = 0
    foldid <- sample(rep(seq(5), length = sum(A == 0)))
    fit_A0 <- cv.glmnet(x = as.matrix(data.table(W[A == 0,])),
                     y = pseudo_outcome_A0, family = "gaussian", weights = pseudo_weights_A0,
                     keep = TRUE, foldid = foldid, alpha = 1, relax = TRUE)
    A0 <- as.numeric(predict(fit_A0, newx = as.matrix(data.table(W)), type = "response"))

    # A = 1
    foldid <- sample(rep(seq(5), length = sum(A == 1)))
    fit_A1 <- cv.glmnet(x = as.matrix(data.table(W[A == 1,])),
                     y = pseudo_outcome_A1, family = "gaussian", weights = pseudo_weights_A1,
                     keep = TRUE, foldid = foldid, alpha = 1, relax = TRUE)
    A1 <- as.numeric(predict(fit_A1, newx = as.matrix(data.table(W)), type = "response"))

    # design matrices
    non_zero_A0 <- which(as.numeric(coef(fit_A0)) != 0)
    non_zero_A1 <- which(as.numeric(coef(fit_A1)) != 0)
    non_zero <- union(non_zero_A0, non_zero_A1)
    x_basis <- cbind(1, as.matrix(data.table(W, A = A)))#[, non_zero, drop = FALSE]
    x_basis_A1 <- cbind(1, as.matrix(data.table(W, A = 1)))#[, non_zero_A1, drop = FALSE]
    x_basis_A0 <- cbind(1, as.matrix(data.table(W, A = 0)))#[, non_zero_A0, drop = FALSE]
  } else if (method == "HAL") {
    fit <- fit_relaxed_hal(as.matrix(data.table(W, A = A)),
                           pseudo_outcome, "gaussian", weights = pseudo_weights)
    pred <- as.vector(fit$pred)

    x_basis <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W, A)))
    x_basis_A1 <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W, A = 1)))
    x_basis_A0 <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W, A = 0)))
    A1 <- as.vector(x_basis_A1 %*% fit$beta)
    A0 <- as.vector(x_basis_A0 %*% fit$beta)
  } else if (method == "glm") {
    # A = 0
    fit_A0 <- lm(pseudo_outcome_A0 ~ ., weights = pseudo_weights_A0, data = data.frame(W[A == 0,]))
    A0 <- as.numeric(predict(fit_A0, newdata = data.frame(W)))

    # A = 1
    fit_A1 <- lm(pseudo_outcome_A1 ~ ., weights = pseudo_weights_A1, data = data.frame(W[A == 1,]))
    A1 <- as.numeric(predict(fit_A1, newdata = data.frame(W)))

    x_basis <- cbind(1, as.matrix(data.table(W, A = A)))
    x_basis_A1 <- cbind(1, as.matrix(data.table(W, A = 1)))
    x_basis_A0 <- cbind(1, as.matrix(data.table(W, A = 0)))
  }

  return(list(A1 = A1,
              A0 = A0,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0))
}

# function to learn psi_tilde, R-learner
learn_psi_tilde <- function(W, A, Y, g, theta, method = "lasso") {
  weights <- 1
  pseudo_outcome <- (Y-theta)/(A-g)
  pseudo_weights <- (A-g)^2*weights

  pred <- NULL
  x_basis <- NULL
  if (method == "lasso") {
    foldid <- sample(rep(seq(5), length = length(A)))
    fit <- cv.glmnet(x = as.matrix(data.table(W)),
                     y = pseudo_outcome, family = "gaussian", weights = pseudo_weights,
                     keep = TRUE, foldid = foldid, alpha = 1, relax = TRUE)
    y_lambda_min <- fit$lambda.min
    pred <- as.numeric(predict(fit, newx = as.matrix(data.table(W)), s = y_lambda_min, type = "response"))

    # design matrix
    x_basis <- cbind(1, as.matrix(data.table(W)))
  } else if (method == "HAL") {
    fit <- fit_relaxed_hal(as.matrix(data.table(W)), pseudo_outcome,
                           "gaussian", weights = pseudo_weights)
    pred <- as.numeric(fit$pred)
    x_basis <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W)))
  } else if (method == "glm") {
    fit <- lm(pseudo_outcome ~ ., data = data.frame(W))
    pred <- as.numeric(predict(fit))
    x_basis <- cbind(1, as.matrix(data.table(W)))
  }

  return(list(pred = pred,
              x_basis = x_basis))
}

# TODO: function to learn psi_tilde using TMLE
learn_psi_tilde_tmle <- function(g, Pi) {
  # clever covariates
  H1_n <- A/g*tau_pred$A1
  H0_n <- (1-A)/(1-g_pred)*tau_pred$A0

  # logistic submodel
  submodel <- glm(S ~ -1 + H1_n + H0_n + offset(Pi_pred$A1), family = "binomial")
  epsilon <- coef(submodel)

  # TMLE updates
  Pi_A1_star <- plogis(Pi_pred$A1 + epsilon["H1_n"] * H1_n)
  Pi_A0_star <- plogis(Pi_pred$A0 + epsilon["H0_n"] * H0_n)
  pred <- A * Pi_A1_star + (1 - A) * Pi_A0_star

  return(list(pred = as.vector(pred),
              A1 = as.vector(Pi_A1_star),
              A0 = as.vector(Pi_A0_star)))
}

# function to learn g(1|S,W)=P(1|S,W)
learn_g <- function(S, W, A, p_rct, method = "lasso") {
  pred <- vector(length = length(A))
  if (method == "lasso") {
    foldid = sample(rep(seq(5), length = sum(S == 0)))
    fit <- cv.glmnet(x = as.matrix(data.table(W[S == 0,])), y = A[S == 0],
                     foldid = foldid, keep = TRUE, alpha = 1, family = "binomial")
    y_lambda_min <- fit$lambda.min
    pred[S == 1] <- rep(p_rct, sum(S))
    pred[S == 0] <- as.numeric(predict(fit, newx = as.matrix(data.table(W[S == 0,])), s = y_lambda_min, type = "response"))
    #pred[S == 0] <- plogis(fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == y_lambda_min])
  } else if (method == "HAL") {
    fit <- fit_hal(X = as.matrix(data.table(W[S == 0,])), Y = A[S == 0],
                   family = "binomial", smoothness_orders = 0)
    pred[S == 1] <- rep(p_rct, sum(S))
    pred[S == 0] <- as.numeric(predict(fit, new_data = as.matrix(data.table(W[S == 0,]))))
  } else if (method == "glm") {
    fit <- glm(A[S == 0] ~ ., data = data.frame(W[S == 0,]), family = "binomial")
    pred[S == 1] <- rep(p_rct, sum(S))
    pred[S == 0] <- as.numeric(predict(fit, type = "response"))
  }

  return(pred)
}

# function to learn g(1|W)=P(1|W)
learn_g_ignore_S <- function(W, A, method = "lasso") {
  pred <- NULL
  if (method == "lasso") {
    foldid = sample(rep(seq(5), length = length(A)))
    fit <- cv.glmnet(x = as.matrix(data.table(W)), y = A,
                     foldid = foldid, keep = TRUE, alpha = 1, family = "binomial")
    y_lambda_min <- fit$lambda.min
    pred <- plogis(fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == y_lambda_min])
  } else if (method == "HAL") {
    fit_g <- fit_relaxed_hal(as.matrix(data.table(W)), A, "binomial")
    pred <- as.vector(fit_g$pred)
  }

  return(pred)
}

learn_g_tmp <- function(W, A, method = "lasso") {
  pred <- NULL
  if (method == "lasso") {
    foldid = sample(rep(seq(5), length = length(A)))
    fit <- cv.glmnet(x = as.matrix(data.frame(W)), y = A,
                     foldid = foldid, keep = TRUE, alpha = 1, family = "binomial")
    y_lambda_min <- fit$lambda.min
    pred <- plogis(fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == y_lambda_min])
  } else if (method == "HAL") {
    fit <- fit_relaxed_hal(as.matrix(data.frame(W)), A, "binomial")
    pred <- as.vector(fit_theta$pred)
  }

  return(pred)
}

# function to perform TMLE update of Pi
Pi_tmle <- function(S, W, A, g, tau, Pi, target_gwt=FALSE) {

  # whether to target using weight
  if(target_gwt){
    wt <- A/g*tau$A1+(1-A)/(1-g)*tau$A0
    H1_n <- A
    H0_n <- 1-A
  } else{
    wt <- rep(1, length(A))
    H1_n <- A/g*tau$A1
    H0_n <- (1-A)/(1-g)*tau$A0
  }

  # logistic submodel
  suppressWarnings(
    epsilon <- coef(glm(S ~ -1 + offset(Pi$pred) + H0_n + H1_n, family = "binomial", weights = wt))
  )
  epsilon[is.na(epsilon)] <- 0

  # TMLE updates
  Pi_star <- NULL
  if (target_gwt) {
    Pi_star$pred <- plogis(Pi$pred+epsilon[1]*H0_n+epsilon[2]*H1_n)
    Pi_star$A0 <- plogis(Pi$A0+epsilon[1])
    Pi_star$A1 <- plogis(Pi$A1+epsilon[2])
  } else {
    Pi_star$pred <- plogis(Pi$pred+epsilon[1]*H0_n+epsilon[2]*H1_n)
    Pi_star$A0 <- plogis(Pi$A0+epsilon[1]/(1-g)*tau$A0)
    Pi_star$A1 <- plogis(Pi$A1+epsilon[2]/g*tau$A1)
  }

  return(Pi_star)
}

get_eic_Pi <- function(g, tau, Pi, S, A) {
  return((A/g*tau$A1-(1-A)/(1-g)*tau$A0)*(S-Pi$pred)-mean(Pi$pred))
}

get_eic_psi_pound_parametric <- function(Pi, tau, g, psi_pound, S, A, Y, n) {
  W_comp <- as.vector((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
  Pi_comp <- as.vector(A/g*tau$A1*(S-Pi$A1)-(1-A)/(1-g)*tau$A0*(S-Pi$A0))
  D_beta <- solve(t(tau$x_basis)%*%tau$x_basis/n)%*%t(tau$x_basis)%*%diag(Y-tau$pred)
  beta_comp <- as.vector((1-Pi$A0)%*%((tau$x_basis_S1A0-tau$x_basis_S0A0)/n)%*%D_beta-
                           (1-Pi$A1)%*%((tau$x_basis_S1A1-tau$x_basis_S0A1)/n)%*%D_beta)
  return(W_comp+Pi_comp+beta_comp-psi_pound)
}

get_eic_psi_pound <- function(Pi, tau, g, theta, psi_pound, S, A, Y, n) {
  W_comp <- (1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1-psi_pound # solved
  Pi_comp <- ((A/g*tau$A1-(1-A)/(1-g)*tau$A0))*(S-Pi$pred) # solved
  IM <- solve(t(tau$x_basis)%*%diag((Pi$pred*(1-Pi$pred)))%*%tau$x_basis/n)
  #IM_A1 <- solve(t(tau$x_basis_A1)%*%diag((Pi$A1*(1-Pi$A1)))%*%tau$x_basis_A1/n)
  #IM_A0 <- solve(t(tau$x_basis_A0)%*%diag((Pi$A0*(1-Pi$A0)))%*%tau$x_basis_A0/n)
  IM_A0 <- IM%*%colMeans(tau$x_basis_A0*(1-Pi$A0))
  IM_A1 <- IM%*%colMeans(tau$x_basis_A1*(1-Pi$A1))
  tmp <- vector(length = n)
  tmp[A == 1] <- tau$A1[A == 1]
  tmp[A == 0] <- tau$A0[A == 0]
  D_beta_A0 <- as.numeric(tau$x_basis%*%IM_A0)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tmp)
  D_beta_A1 <- as.numeric(tau$x_basis%*%IM_A1)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tmp)
  beta_comp <- D_beta_A0-D_beta_A1
  #D_beta_A1 <- as.numeric(tau$x_basis_A1%*%IM_A1)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tmp)
  #D_beta_A0 <- as.numeric(tau$x_basis_A0%*%IM_A0)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tmp)
  #beta_comp <- D_beta_A0*(1-Pi$A0)-D_beta_A1*(1-Pi$A1)

  return(W_comp+Pi_comp+beta_comp)
}

get_eic_psi_tilde <- function(psi_tilde, g_pred, theta, Y, A, n) {
  IM <- solve(t(psi_tilde$x_basis)%*%diag((g_pred*(1-g_pred)))%*%psi_tilde$x_basis/n)%*%colMeans(psi_tilde$x_basis)
  D_beta <- psi_tilde$x_basis%*%IM*(A-g_pred)*(Y-theta-(A-g_pred)*psi_tilde$pred)
  return(as.vector(psi_tilde$pred-mean(psi_tilde$pred)+D_beta))
}

learn_Q <- function(W, A, Y) {
  fit_Q <- fit_relaxed_hal(as.matrix(data.table(W, A)), Y, "gaussian")
  pred <- as.vector(fit_Q$pred)

  x_basis_A1 <- make_counter_design_matrix(fit_Q$basis_list, as.matrix(data.table(W, A = 1)))
  x_basis_A0 <- make_counter_design_matrix(fit_Q$basis_list, as.matrix(data.table(W, A = 0)))
  A1 <- as.vector(x_basis_A1 %*% fit_Q$beta)
  A0 <- as.vector(x_basis_A0 %*% fit_Q$beta)

  return(list(pred = pred,
              A1 = A1,
              A0 = A0))
}

Q_tmle <- function(g, Q, A, Y_bound) {

  wt <- A/g+(1-A)/(1-g)
  H1W <- A
  H0W <- 1-A

  submodel <- glm(Y_bound ~ -1 + offset(Q$pred) + H0W + H1W, family = "quasibinomial", weights = wt)
  epsilon <- coef(submodel)

  Q_star <- Q$pred+epsilon[1]*H0W+epsilon[2]*H1W
  Q_A1_star <- Q$A1+rep(epsilon[1], length(Y_bound))
  Q_A0_star <- Q$A0+rep(epsilon[2], length(Y_bound))

  return(list(Q_star = Q_star,
              A1 = Q_A1_star,
              A0 = Q_A0_star))
}

get_eic_psi_tilde_2 <- function(g, Q_A1, Q_A0, A, Y) {
  D_A1 <- A/g*(Y-Q_A1)
  D_A0 <- (1-A)/(1-g)*(Y-Q_A0)

  return(D_A1-D_A0+Q_A1-Q_A0-mean(Q_A1-Q_A0))
}

bound <- function(X) {
  X_max <- max(X, na.rm = TRUE)
  X_min <- min(X, na.rm = TRUE)

  return((X-X_min)/(X_max-X_min))
}

learn_tau_parametric <- function(S, W, A, Y) {
  fit_tau <- fit_relaxed_hal(X = as.matrix(data.table(S, W, A)), Y = Y, family = "gaussian",
                             smoothness_orders = 1, num_knots = c(50,50,20))
  pred <- as.vector(fit_tau$pred)

  x_basis <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S, W, A)))
  x_basis_A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S, W, A = 1)))
  x_basis_A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S, W, A = 0)))
  x_basis_S1A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S = 1, W, A = 1)))
  x_basis_S0A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S = 0, W, A = 1)))
  x_basis_S1A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S = 1, W, A = 0)))
  x_basis_S0A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S = 0, W, A = 0)))

  A1 <- as.vector(x_basis_S1A1 %*% fit_tau$beta) - as.vector(x_basis_S0A1 %*% fit_tau$beta)
  A0 <- as.vector(x_basis_S1A0 %*% fit_tau$beta) - as.vector(x_basis_S0A0 %*% fit_tau$beta)

  return(list(pred = pred,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0,
              x_basis_S1A1 = x_basis_S1A1,
              x_basis_S0A1 = x_basis_S0A1,
              x_basis_S1A0 = x_basis_S1A0,
              x_basis_S0A0 = x_basis_S0A0,
              A1 = A1,
              A0 = A0))
}

# psi_pound_pred <- mean((1-Pi_star$Pi_A0_star)*tau_star$tau_A0_star)-mean((1-Pi_star$Pi_A1_star)*tau_star$tau_A1_star)
#

#
# # Psi_tilde --------------------------------------------------------------------
#
# # learn theta_tilde(W)=E(Y|W)
# fit_theta_tilde <- fit_relaxed_hal(as.matrix(data.table(W)), Y, "gaussian")
# theta_tilde_pred <- fit_theta_tilde$pred
#
# # learn g(1|W)=P(A=1|W)
# fit_g <- fit_relaxed_hal(as.matrix(data.table(W)), A, "binomial")
# g_pred <- fit_g$pred
#
# weights <- 1
#
# # make pseudo outcome and weights
# pseudo_outcome_tilde <- ifelse(abs(A - g_pred) < 1e-10, 0, (Y - theta_tilde_pred) / (A - g_pred))
# pseudo_weights_tilde <- (A - g_pred)^2 * weights
#
# # learn Psi_tilde
# keep <- which(abs(A - g_pred) > 1e-10)
# fit_psi_tilde <- fit_relaxed_hal(as.matrix(data.table(W)), pseudo_outcome_tilde[keep], "gaussian", weights = pseudo_weights_tilde[keep])
# psi_tilde_pred <- mean(fit_psi_tilde$pred)
#
# # target parameter point estimate
# x_basis_A1 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 1)))
# x_basis_A0 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 0)))
# Pi_A1_pred <- x_basis_A1 %*% fit_Pi$beta
# Pi_A0_pred <- x_basis_A0 %*% fit_Pi$beta
# psi_pound_pred <- mean((1 - Pi_A0_pred) * tau_A0_pred - (1 - Pi_A1_pred) * tau_A1_pred)
# psi_pred <- psi_tilde_pred - psi_pound_pred
# print(psi_pred)
