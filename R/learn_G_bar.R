learn_G_bar <- function(data,
                        W,
                        A,
                        g,
                        method,
                        folds,
                        G_bar_bounds) {

  if (method == "glm") {
    data_A1 <- copy(data); data_A1[, A := 1]
    data_A0 <- copy(data); data_A0[, A := 0]

    fit <- glm(data[["Delta_G"]] ~ ., data = data, family = "binomial")
    A1 <- predict(fit, newdata = data_A1, type = "response")
    A0 <- predict(fit, newdata = data_A0, type = "response")
  }

  return(list(A1 = .bound(as.numeric(A1), G_bar_bounds),
              A0 = .bound(as.numeric(A0), G_bar_bounds),
              integrate_A = as.numeric(A1)*g+as.numeric(A0)*(1-g)))
}
