Q_tmle <- function(g, Q, A, Y_bound) {
  wt <- A / g + (1 - A) / (1 - g)
  H1W <- A
  H0W <- 1 - A

  submodel <- glm(Y_bound ~ -1 + offset(Q$pred) + H0W + H1W, family = "quasibinomial", weights = wt)
  epsilon <- coef(submodel)

  Q_star <- Q$pred + epsilon[1] * H0W + epsilon[2] * H1W
  Q_A1_star <- Q$A1 + rep(epsilon[1], length(Y_bound))
  Q_A0_star <- Q$A0 + rep(epsilon[2], length(Y_bound))

  return(list(
    Q_star = Q_star,
    A1 = Q_A1_star,
    A0 = Q_A0_star
  ))
}

.bound <- function(x, bounds) {
  return(pmin(pmax(x, bounds[1]), bounds[2]))
}

#' @import sl3
learn_SW <- function(S,
                     W,
                     folds,
                     method,
                     Pi_bounds) {

  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners(family)
  }

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                         Lrnr_cv_selector$new(loss_loglik_binomial))
    task_SW <- sl3_Task$new(data = data.table(S=S, W),
                            covariates = colnames(W),
                            outcome = "S",
                            folds = folds,
                            outcome_type = "binomial")
    suppressMessages(fit_SW <- lrnr$train(task_SW))
    pred <- fit_SW$predict(task_SW)
  } else if (method == "glm") {
    fit <- glm(S ~ ., data = data.frame(S=S, W), family = "binomial")
    pred <- as.numeric(predict(fit, newdata = data.frame(W), type = "response"))
  }

  return(.bound(pred, Pi_bounds))
}

target_Q <- function(S,
                     W,
                     A,
                     Y,
                     Pi,
                     g11W,
                     Q) {
  # bound Y
  min_Y <- min(Y, Q$Q1WA, Q$Q1W1, Q$Q1W0)-0.001
  max_Y <- max(Y, Q$Q1WA, Q$Q1W1, Q$Q1W0)+0.001
  Y_scaled <- (Y-min_Y)/(max_Y-min_Y)
  Q1WA_scaled <- (Q$Q1WA-min_Y)/(max_Y-min_Y)
  Q1W1_scaled <- (Q$Q1W1-min_Y)/(max_Y-min_Y)
  Q1W0_scaled <- (Q$Q1W0-min_Y)/(max_Y-min_Y)

  # clever covariate
  HSAW <- (S/Pi)*(A/g11W-(1-A)/(1-g11W))
  HS1W <- (S/Pi)*(1/g11W)
  HS0W <- (S/Pi)*(-1/(1-g11W))

  # logistic submodel
  epsilon <- coef(glm(Y_scaled ~ -1+offset(qlogis(Q1WA_scaled))+HSAW,
                      family = "quasibinomial"))
  epsilon[is.na(epsilon)] <- 0

  # TMLE updates
  Q_star <- list(Q1WA = plogis(qlogis(Q1WA_scaled)+epsilon*HSAW),
                 Q1W1 = plogis(qlogis(Q1W1_scaled)+epsilon*HS1W),
                 Q1W0 = plogis(qlogis(Q1W0_scaled)+epsilon*HS0W))

  # scale back
  Q_star$Q1WA <- Q_star$Q1WA*(max_Y-min_Y)+min_Y
  Q_star$Q1W1 <- Q_star$Q1W1*(max_Y-min_Y)+min_Y
  Q_star$Q1W0 <- Q_star$Q1W0*(max_Y-min_Y)+min_Y

  return(Q_star)
}

target_Q_rct_W <- function(S,
                           W,
                           A,
                           Y,
                           pS,
                           g11W,
                           Q) {
  # bound Y
  min_Y <- min(Y, Q$Q1WA, Q$Q1W1, Q$Q1W0)-0.001
  max_Y <- max(Y, Q$Q1WA, Q$Q1W1, Q$Q1W0)+0.001
  Y_scaled <- (Y-min_Y)/(max_Y-min_Y)
  Q1WA_scaled <- (Q$Q1WA-min_Y)/(max_Y-min_Y)
  Q1W1_scaled <- (Q$Q1W1-min_Y)/(max_Y-min_Y)
  Q1W0_scaled <- (Q$Q1W0-min_Y)/(max_Y-min_Y)

  # clever covariate
  HSAW <- (S/pS)*(A/g11W-(1-A)/(1-g11W))
  HS1W <- (S/pS)*(1/g11W)
  HS0W <- (S/pS)*(-1/(1-g11W))

  # logistic submodel
  epsilon <- coef(glm(Y_scaled ~ -1+offset(qlogis(Q1WA_scaled))+HSAW,
                      family = "quasibinomial"))
  epsilon[is.na(epsilon)] <- 0

  # TMLE updates
  Q_star <- list(Q1WA = plogis(qlogis(Q1WA_scaled)+epsilon*HSAW),
                 Q1W1 = plogis(qlogis(Q1W1_scaled)+epsilon*HS1W),
                 Q1W0 = plogis(qlogis(Q1W0_scaled)+epsilon*HS0W))

  # scale back
  Q_star$Q1WA <- Q_star$Q1WA*(max_Y-min_Y)+min_Y
  Q_star$Q1W1 <- Q_star$Q1W1*(max_Y-min_Y)+min_Y
  Q_star$Q1W0 <- Q_star$Q1W0*(max_Y-min_Y)+min_Y

  return(Q_star)
}

#' @import sl3
#' @import data.table
learn_QSWA <- function(S,
                       W,
                       A,
                       Y,
                       folds,
                       folds_S1,
                       family,
                       method,
                       pooling) {

  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners(family)
  }

  if (is.list(method)) {
    if (family == "binomial") {
      loss <- loss_loglik_binomial
    } else if (family == "gaussian") {
      loss <- loss_squared_error
    }
    lrnr_stack <- Stack$new(method)
    lrnr <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                         Lrnr_cv_selector$new(loss))
    if (pooling) {
      task_QSWA <- sl3_Task$new(data = data.table(Y=Y, S=S, A=A, W),
                                covariates = c(colnames(W), "S", "A"),
                                outcome = "Y",
                                folds = folds,
                                outcome_type = family)
      task_Q1WA <- sl3_Task$new(data = data.table(Y=Y, S=1, A=A, W),
                                covariates = c(colnames(W), "S", "A"),
                                outcome = "Y",
                                folds = folds,
                                outcome_type = family)
      task_Q1W1 <- sl3_Task$new(data = data.table(Y=Y, S=1, A=1, W),
                                covariates = c(colnames(W), "S", "A"),
                                outcome = "Y",
                                folds = folds,
                                outcome_type = family)
      task_Q1W0 <- sl3_Task$new(data = data.table(Y=Y, S=1, A=0, W),
                                covariates = c(colnames(W), "S", "A"),
                                outcome = "Y",
                                folds = folds,
                                outcome_type = family)
    } else {
      task_QSWA <- sl3_Task$new(data = data.table(Y=Y, A=A, W)[S == 1,],
                                covariates = c(colnames(W), "A"),
                                outcome = "Y",
                                folds = folds_S1,
                                outcome_type = family)
      task_Q1WA <- sl3_Task$new(data = data.table(Y=Y, A=A, W),
                                covariates = c(colnames(W), "A"),
                                outcome = "Y",
                                folds = folds,
                                outcome_type = family)
      task_Q1W1 <- sl3_Task$new(data = data.table(Y=Y, A=1, W),
                                covariates = c(colnames(W), "A"),
                                outcome = "Y",
                                folds = folds,
                                outcome_type = family)
      task_Q1W0 <- sl3_Task$new(data = data.table(Y=Y, A=0, W),
                                covariates = c(colnames(W), "A"),
                                outcome = "Y",
                                folds = folds,
                                outcome_type = family)
    }
    suppressMessages(fit_QSWA <- lrnr$train(task_QSWA))
    pred_Q1WA <- fit_QSWA$predict(task_Q1WA)
    pred_Q1W1 <- fit_QSWA$predict(task_Q1W1)
    pred_Q1W0 <- fit_QSWA$predict(task_Q1W0)
  } else if (method == "glm") {
    if (pooling) {
      fit <- glm(Y ~ .^3, data = data.frame(S=S, W, A=A, Y=Y), family = family)
      pred_Q1WA <- as.numeric(predict(fit, newdata = data.frame(S=1, W, A=A), type = "response"))
      pred_Q1W1 <- as.numeric(predict(fit, newdata = data.frame(S=1, W, A=1), type = "response"))
      pred_Q1W0 <- as.numeric(predict(fit, newdata = data.frame(S=1, W, A=0), type = "response"))
    } else {
      fit <- glm(Y ~ .^2, data = data.frame(W, A=A, Y=Y)[S == 1,], family = family)
      pred_Q1WA <- as.numeric(predict(fit, newdata = data.frame(W, A=A, Y=Y), type = "response"))
      pred_Q1W1 <- as.numeric(predict(fit, newdata = data.frame(W, A=1, Y=Y), type = "response"))
      pred_Q1W0 <- as.numeric(predict(fit, newdata = data.frame(W, A=0, Y=Y), type = "response"))
    }
  }

  return(list(Q1WA = pred_Q1WA,
              Q1W1 = pred_Q1W1,
              Q1W0 = pred_Q1W0))
}

learn_g11W <- function(S,
                       W,
                       A,
                       method,
                       g_bounds) {

  if (method == "glm") {
    fit <- glm(A[S == 1] ~ ., data = W[S == 1, , drop = FALSE], family = "binomial")
    pred <- as.numeric(predict(fit, newdata = W, type = "response"))
  }

  return(.bound(pred, g_bounds))
}
