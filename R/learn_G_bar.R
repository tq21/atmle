learn_G_bar <- function(W,
                        A,
                        Delta_G,
                        g,
                        method,
                        folds,
                        G_bar_bounds) {

  if (method == "glm") {
    X <- cbind(W, A = A)
    X_A1 <- cbind(W, A = 1)
    X_A0 <- cbind(W, A = 0)

    fit <- glm(Delta_G ~ ., data = X, family = "binomial")
    A1 <- as.numeric(predict(fit, newdata = X_A1, type = "response"))
    A0 <- as.numeric(predict(fit, newdata = X_A0, type = "response"))
  }

  return(list(A1 = A1,
              A0 = A0,
              integrate_A = A1*g+A0*(1-g)))
}
