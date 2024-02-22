#' @title Fit relaxed HAL
#'
#' @param X A numeric matrix of covariates.
#' @param Y A numeric vector of outcomes.
#' @param family A character string specifying the family of the outcome.
#' @param weights A numeric vector of weights.
#' @param X_no_basis A numeric matrix of covariates that should not be used
#' to construct HAL basis.
#' @param X_unpenalized A numeric matrix of covariates that should not be
#' penalized by HAL.
#'
#' @return A list with the following components:
#' \item{beta}{A numeric vector of coefficients.}
#' \item{pred}{A numeric vector of predictions.}
#' \item{hal_basis_list}{A list of HAL basis.}
#' \item{x_basis}{A numeric matrix of basis.}
fit_relaxed_hal <- function(X,
                            Y,
                            X_unpenalized,
                            X_weak_penalized,
                            X_weak_penalized_level,
                            family,
                            weights,
                            relaxed,
                            v_folds,
                            hal_args) {

  x_basis <- NULL
  hal_basis_list <- NULL

  # use lasso to screen main terms and interactions
  # then unpenalize them in HAL fit
  # if (screen_unpenalize) {
  #   # make design matrix for main terms and interactions
  #   max_degree <- ifelse(ncol(X) >= 20, 2, 3)
  #   X_unpenalized <- model.matrix(as.formula(paste0("Y ~ -1 + (.)^", max_degree)),
  #                                 data = data.frame(X))
  #   screen_fit <- cv.glmnet(x = X_unpenalized, y = Y, family = family,
  #                           weights = weights, alpha = 1, nfolds = v_folds)
#
  #   # get indices of selected terms, excluding intercept
  #   non_zero <- which(as.numeric(coef(screen_fit, s = "lambda.1se")) != 0)[-1] - 1
#
  #   if (length(non_zero) == 0) {
  #     X_unpenalized <- NULL
  #   } else {
  #     X_unpenalized <- X_unpenalized[, non_zero, drop = FALSE]
  #   }
  # }

  # make basis list
  basis_list <- enumerate_basis(x = X,
                                max_degree = ifelse(ncol(X) >= 20, 2, 3),
                                smoothness_orders = 1,
                                num_knots = hal9001:::num_knots_generator(
                                  max_degree = ifelse(ncol(X) >= 20, 2, 3),
                                  smoothness_orders = 1,
                                  base_num_knots_0 = 200,
                                  base_num_knots_1 = 50))
  penalty_factor <- rep(1, length(basis_list))

  if (!is.null(X_weak_penalized)) {
    # make basis list for weakly penalized terms
    # append to basis list
    X_min <- apply(X_weak_penalized, 2, min)
    basis_list_main_terms <- enumerate_basis(x = X_min,
                                             max_degree = ifelse(ncol(X) >= 20, 2, 3),
                                             smoothness_orders = 1)
    basis_list <- c(basis_list, basis_list_main_terms)
    penalty_factor <- c(penalty_factor, rep(X_weak_penalized_level,
                                            length(basis_list_main_terms)))
  }

  # default set of arguments for HAL
  hal_args_default <- list(X = X,
                           Y = Y,
                           formula = NULL,
                           X_unpenalized = X_unpenalized,
                           basis_list = basis_list,
                           penalty_factor = penalty_factor,
                           reduce_basis = NULL,
                           family = family,
                           weights = weights,
                           #fit_control = list(use_min = FALSE),
                           return_x_basis = TRUE)

  # merge arguments
  hal_args <- modifyList(hal_args_default, hal_args)

  # fit HAL
  hal_fit <- do.call(fit_hal, hal_args)
  hal_basis_list <- list()

  if (!is.null(X_unpenalized)) {
    hal_bases_idx <- 2:(length(hal_fit$coefs)-ncol(X_unpenalized))
    hal_basis_list <- hal_fit$basis_list[hal_fit$coefs[hal_bases_idx] != 0]
  } else {
    hal_basis_list <- hal_fit$basis_list[hal_fit$coefs[-1] != 0]
  }

  selected_idx <- which(hal_fit$coefs != 0)
  selected_unpenalized_idx <- which(hal_fit$coefs[length(hal_fit$basis_list)+1:length(hal_fit$coefs)] != 0)-1
  x_basis <- cbind(1, hal_fit$x_basis)[, selected_idx, drop = FALSE]
  beta <- NULL

  if (relaxed) {
    # relaxed HAL fit
    if (family == "binomial") {
      hal_relaxed_fit <- glm.fit(x = x_basis, y = Y, family = binomial(), weights = weights, intercept = FALSE)
      beta <- coef(hal_relaxed_fit)
      beta[is.na(beta)] <- 0
      pred <- x_basis %*% beta
      pred <- as.vector(1 / (1 + exp(-pred)))
    } else if (family == "gaussian") {
      hal_relaxed_fit <- glm.fit(x = x_basis, y = Y, family = gaussian(), weights = weights, intercept = FALSE)
      beta <- coef(hal_relaxed_fit)
      beta[is.na(beta)] <- 0
      pred <- as.vector(x_basis %*% beta)
    }
  } else {
    # regular HAL fit
    beta <- hal_fit$coefs[selected_idx]
    pred <- predict(hal_fit,
                    new_data = X,
                    new_X_unpenalized = X_unpenalized,
                    type = "response")
  }

  #print("number of non-zero coefficients: " %+% length(beta))

  return(list(beta = beta,
              hal_basis_list = hal_basis_list,
              pred = pred,
              x_basis = as.matrix(x_basis),
              selected_unpenalized_idx = selected_unpenalized_idx,
              hal_fit = hal_fit))
}
