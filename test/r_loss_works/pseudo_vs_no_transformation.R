source("utils.R")
learn_tau_test <- function(S, W, A, Y, Pi, theta, method = "lasso", transform=TRUE) {
  pred <- numeric(length = length(A))
  A1 <- numeric(length = length(A))
  A0 <- numeric(length = length(A))
  x_basis <- NULL
  x_basis_A1 <- NULL
  x_basis_A0 <- NULL

  X <- (S-Pi$pred) * data.frame(W, A)
  X_A1 <- (S-Pi$pred) * data.frame(W, A = 1)
  X_A0 <- (S-Pi$pred) * data.frame(W, A = 0)

  if (method == "lasso") {
    if (transform) {
      set.seed(123)
      weights <- 1
      pseudo_outcome <- (Y-theta)/(S-Pi$pred)
      pseudo_weights <- (S-Pi$pred)^2*weights
      fit <- cv.glmnet(x = as.matrix(X), y = pseudo_outcome, weights = pseudo_weights,
                       family = "gaussian", keep = TRUE, nfolds = 5, alpha = 1, relax = TRUE)
      pred <- as.numeric(predict(fit, newx = as.matrix(X), type = "response"))
      A1 <- as.numeric(predict(fit, newx = as.matrix(X_A1), type = "response"))
      A0 <- as.numeric(predict(fit, newx = as.matrix(X_A0), type = "response"))
    } else {
      set.seed(123)
      pseudo_outcome <- Y-theta
      fit <- cv.glmnet(x = as.matrix(X), y = pseudo_outcome,
                       family = "gaussian", keep = TRUE, nfolds = 5, alpha = 1, relax = TRUE)
      pred <- as.numeric(predict(fit, newx = as.matrix(X), type = "response"))
      A1 <- as.numeric(predict(fit, newx = as.matrix(X_A1), type = "response"))
      A0 <- as.numeric(predict(fit, newx = as.matrix(X_A0), type = "response"))

      # design matrices
      non_zero <- which(as.numeric(coef(fit)) != 0)
      x_basis <- cbind(1, as.matrix(data.table(W, A = A)))[, non_zero, drop = FALSE]
      x_basis_A1 <- cbind(1, as.matrix(data.table(W, A = 1)))[, non_zero, drop = FALSE]
      x_basis_A0 <- cbind(1, as.matrix(data.table(W, A = 0)))[, non_zero, drop = FALSE]
    }

  } else if (method == "HAL") {
    pseudo_outcome <- (Y-theta)/(S-Pi$pred)
    pseudo_weights <- (S-Pi$pred)^2*weights
    fit <- fit_relaxed_hal(as.matrix(data.table(W, A = A)),
                           pseudo_outcome, "gaussian", weights = pseudo_weights)
    pred <- as.vector(fit$pred)

    x_basis <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W, A)))
    x_basis_A1 <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W, A = 1)))
    x_basis_A0 <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W, A = 0)))
    A1 <- as.vector(x_basis_A1 %*% fit$beta)
    A0 <- as.vector(x_basis_A0 %*% fit$beta)
  } else if (method == "glm") {
    # A = 1
    fit_A1 <- lm(pseudo_outcome_A1 ~ ., weights = pseudo_weights_A1, data = X_A1)
    A1 <- as.numeric(predict(fit_A1, newdata = X))

    # A = 0
    fit_A0 <- lm(pseudo_outcome_A0 ~ ., weights = pseudo_weights_A0, data = X_A0)
    A0 <- as.numeric(predict(fit_A0, newdata = X))

    x_basis <- cbind(1, as.matrix(data.table(W, A = A)))
    x_basis_A1 <- cbind(1, as.matrix(data.table(W, A = 1)))
    x_basis_A0 <- cbind(1, as.matrix(data.table(W, A = 0)))
  }

  return(list(A0 = A0))
}


B <- 200
res_transform <- rep(NA, B)
res_no_transform <- rep(NA, B)
for (i in 1:B) {
  print(i)
  data <- generate_data(500, 1.5, 0)
  S_node <- 1
  W_node <- c(2, 3, 4, 5)
  A_node <- 6
  Y_node <- 7
  nuisance_method <- "lasso"
  working_model <- "lasso"
  p_rct=0.5
  verbose=TRUE

  S <- data[, S_node]
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]
  n <- nrow(data)
  theta <- learn_theta(W, A, Y, nuisance_method)
  Pi <- learn_Pi(S, W, A, nuisance_method)
  g <- learn_g(S, W, A, p_rct, nuisance_method)
  tau_transform <- learn_tau_test(S, W, A, Y, Pi, theta, method = working_model, transform=TRUE)
  tau_no_transform <- learn_tau_test(S, W, A, Y, Pi, theta, method = working_model, transform=FALSE)
  res_transform[i] <- mean(tau_transform$A0)
  res_no_transform[i] <- mean(tau_no_transform$A0)
}

summary(res_transform)
hist(res_no_transform)
