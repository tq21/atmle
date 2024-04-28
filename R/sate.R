# The following code is adapted from code for the paper
# "targeted estimation and inference for the sample average treatment effect
# in trials with and without pair-matching"
# by Balzer, Petersen, and van der Laan (2016)
sate <- function(data,
                 S_node,
                 W_node,
                 A_node,
                 Y_node,
                 g_rct,
                 family) {
  # define nodes ---------------------------------------------------------------
  S <- data[, S_node] # study indicator
  data <- data[S == 1, ]
  W <- data[, W_node] # covariates
  A <- data[, A_node] # treatment
  Y <- data[, Y_node] # outcome
  n <- nrow(data) # sample size

  # sample average counterfactual outcome under A=a
  EY_1 <- mean(Y[A == 1])
  EY_0 <- mean(Y[A == 0])

  # sample average treatment effect
  sample_ate <- EY_1 - EY_0

  # TMLE for SATE
  res <- doTMLE(W, A, Y, g_rct, family)

  return(res)
}

doTMLE <- function(W, A, Y,
                   g_rct,
                   family) {
  n <- length(A)

  # make counterfactual data
  X <- cbind(W, A = A, Y = Y)
  X1 <- copy(X)
  X0 <- copy(X)
  X1$A <- 1
  X0$A <- 0

  glm.out <- glm(Y ~ ., family = family, data = X)

  # get predicted outcomes under obs exp, txt and control
  QbarAW <- as.numeric(predict(glm.out, newdata = X, type = "response"))
  Qbar1W <- as.numeric(predict(glm.out, newdata = X1, type = "response"))
  Qbar0W <- as.numeric(predict(glm.out, newdata = X0, type = "response"))

  # propensity score
  pscore <- rep(g_rct, n)

  # clever covariates
  H.1W <- 1 / pscore
  H.0W <- -1 / (1 - pscore)
  H.AW <- rep(NA, n)
  H.AW[A == 1] <- H.1W[A == 1]
  H.AW[A == 0] <- H.0W[A == 0]

  # bound Y
  min_Y <- min(Y, QbarAW, Qbar1W, Qbar0W) - 0.001
  max_Y <- max(Y, QbarAW, Qbar1W, Qbar0W) + 0.001
  Y_bounded <- (Y - min_Y) / (max_Y - min_Y)
  QbarAW <- (QbarAW - min_Y) / (max_Y - min_Y)
  Qbar1W <- (Qbar1W - min_Y) / (max_Y - min_Y)
  Qbar0W <- (Qbar0W - min_Y) / (max_Y - min_Y)

  # updating step
  logitUpdate <- suppressWarnings(
    glm(Y_bounded ~ -1 + offset(qlogis(QbarAW)) + H.AW, family = "quasibinomial")
  )

  # estimated coefficient on the clever covariate
  eps <- logitUpdate$coef

  # targeted estimates of the outcome regression
  QbarAW <- plogis(qlogis(QbarAW) + eps * H.AW)
  Qbar0W <- plogis(qlogis(Qbar0W) + eps * H.0W)
  Qbar1W <- plogis(qlogis(Qbar1W) + eps * H.1W)

  # scale back
  QbarAW <- QbarAW * (max_Y - min_Y) + min_Y
  Qbar0W <- Qbar0W * (max_Y - min_Y) + min_Y
  Qbar1W <- Qbar1W * (max_Y - min_Y) + min_Y

  # risk estimates under txt, under control and risk difference
  EY_1 <- mean(Qbar1W)
  EY_0 <- mean(Qbar0W)
  sample_ate <- mean(Qbar1W - Qbar0W)

  # the relevant components of the influence curve
  DY <- H.AW * (Y - QbarAW)
  var.SATE <- var(DY) / n
  lower <- sample_ate - 1.96 * sqrt(var.SATE)
  upper <- sample_ate + 1.96 * sqrt(var.SATE)

  return(data.frame(
    est = sample_ate,
    lower = lower,
    upper = upper
  ))
}
