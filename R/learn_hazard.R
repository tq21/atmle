#' @title Estimate the conditional hazard function for event (failure/censoring)
#'
#' @param data_long A `data.table` of long format data.
#' @param X A `character` vector of covariate names.
#' @param T_tilde A vector of observed failure/censored times.
#' @param event A `character` of event variable name.
#' @param method A `character` of method to use for hazard estimation.
#' @param folds A list of folds for cross-fitting.
#' @param cross_fit_nuisance A `logical` of whether to cross-fit.
#' @param counter_var A `character` of counterfactual variable name, if any.
#'
#' @export
learn_hazard <- function(data_long,
                         X,
                         T_tilde,
                         event,
                         method,
                         folds,
                         cross_fit_nuisance,
                         counter_var = NULL) {

  pred <- A0 <- A1 <- numeric(length = data_long[, .N])

  if (!is.null(counter_var)) {
    # make counterfactual prediction data
    data_long_A1 <- copy(data_long)
    data_long_A1 <- data_long_A1[, (counter_var) := 1]
    data_long_A0 <- copy(data_long)
    data_long_A0 <- data_long_A0[, (counter_var) := 0]
  }

  # make training data
  data_train <- copy(data_long)
  data_train <- data_train[(T_tilde) >= t]
  cov_names <- c(X, "t")

  if (method == "glm") {
    # fit hazard regression using glm
    if (cross_fit_nuisance) {
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(data_train[id %in% train_idx][[event]] ~ .,
                   data = data_train[id %in% train_idx, ..cov_names],
                   family = "binomial")
        pred[which(data_long[["id"]] %in% valid_idx)] <<-
          predict(fit, newdata = data_long[id %in% valid_idx, ..cov_names], type = "response")

        if (!is.null(counter_var)) {
          A1[which(data_long[["id"]] %in% valid_idx)] <<-
            predict(fit, newdata = data_long_A1[id %in% valid_idx, ..cov_names], type = "response")
          A0[which(data_long[["id"]] %in% valid_idx)] <<-
            predict(fit, newdata = data_long_A0[id %in% valid_idx, ..cov_names], type = "response")
        }
      })
    } else {
      fit <- glm(data_train[[event]] ~ .,
                 data = data_train[, ..cov_names],
                 family = "binomial")
      pred <- predict(fit, newdata = data_long[, ..cov_names], type = "response")

      if (!is.null(counter_var)) {
        A1 <- predict(fit, newdata = data_long_A1[, ..cov_names], type = "response")
        A0 <- predict(fit, newdata = data_long_A0[, ..cov_names], type = "response")
      }
    }
  } else if (method == "glmnet") {
    # fit hazard regression using lasso
    if (cross_fit_nuisance) {
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- cv.glmnet(x = as.matrix(data_train[id %in% train_idx, ..cov_names]),
                         y = data_train[id %in% train_idx][[event]],
                         family = "binomial", keep = TRUE, alpha = 1,
                         nfolds = length(folds))
        pred[which(data_long[["id"]] %in% valid_idx)] <<-
          predict(fit, newx = as.matrix(data_long[id %in% valid_idx, ..cov_names]),
                  s = "lambda.min", type = "response")

        if (!is.null(counter_var)) {
          A1[which(data_long[["id"]] %in% valid_idx)] <<-
            predict(fit, newx = as.matrix(data_long_A1[id %in% valid_idx, ..cov_names]),
                    s = "lambda.min", type = "response")
          A0[which(data_long[["id"]] %in% valid_idx)] <<-
            predict(fit, newx = as.matrix(data_long_A0[id %in% valid_idx, ..cov_names]),
                    s = "lambda.min", type = "response")
        }
      })
    } else {
      fit <- cv.glmnet(x = as,matrix(data_train[, ..cov_names]),
                       y = data_train[[event]],
                       family = "binomial", keep = TRUE, alpha = 1,
                       nfolds = length(folds))
      pred <- predict(fit, newx = as.matrix(data_long[, ..cov_names]),
                      s = "lambda.min", type = "response")

      if (!is.null(counter_var)) {
        A1 <- predict(fit, newx = as.matrix(data_long_A1[, ..cov_names]),
                      s = "lambda.min", type = "response")
        A0 <- predict(fit, newx = as.matrix(data_long_A0[, ..cov_names]),
                      s = "lambda.min", type = "response")
      }
    }

  } else if (method == "HAL") {
    # fit hazard regression using HAL
    # currently don't do cross fit here, needs a bit of work to extract
    # fits for each fold and do out-of-sample predictions so to avoid refitting for each fold
    fit <- fit_hal(X = data_train[, ..cov_names],
                   Y = data_train[[event]],
                   id = data_train[["id"]],
                   smoothness_orders = 1,
                   max_degree = 3,
                   family = "binomial",
                   return_x_basis = FALSE)

    # predict hazard
    pred <- predict(fit, new_data = data_long[, ..cov_names], type = "response")

    if (!is.null(counter_var)) {
      A1 <- predict(fit, new_data = data_long_A1[, ..cov_names], type = "response")
      A0 <- predict(fit, new_data = data_long_A0[, ..cov_names], type = "response")
    }
  }

  if (is.null(counter_var)) {
    return(as.numeric(pred))
  } else {
    return(list(pred = as.numeric(pred),
                A1 = as.numeric(A1),
                A0 = as.numeric(A0)))
  }
}
