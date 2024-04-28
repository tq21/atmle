#' @title Relaxed highly adaptive lasso
#'
#' @description First fit highly adaptive lasso, then fit a `glm` only on the
#' design matrix of the selected basis functions (with non-zero coefficients).
#'
#' @param X A numeric matrix of covariates.
#' @param Y A numeric vector of outcomes.
#' @param X_unpenalized A numeric matrix of covariates that should not be
#' penalized by HAL. Used for enforcing minimal working models.
#' @param family A character string specifying the family of the outcome.
#' Either "gaussian" or "binomial".
#' @param weights A numeric vector of weights.
#' @param enumerate_basis_args A list of arguments to pass to `enumerate_basis`
#' function in `hal9001`.
#' @param fit_hal_args A list of arguments to pass to `fit_hal` function in
#' `hal9001`.
#'
#' @importFrom hal9001 enumerate_basis
#' @importFrom hal9001 fit_hal
#'
#' @return A list with the following components:
#' \item{beta}{A numeric vector of coefficients.}
#' \item{hal_basis_list}{A list of selected HAL basis.}
#' \item{x_basis}{The design matrix of selected basis.}
#' \item{pred}{A numeric vector of predictions.}
#' \item{hal_fit}{The output of `fit_hal`.}
fit_relaxed_hal <- function(X,
                            Y,
                            X_unpenalized,
                            X_weak_penalized,
                            X_weak_penalized_level,
                            family,
                            weights,
                            relaxed,
                            enumerate_basis_args = list(),
                            fit_hal_args = list()) {

  # make basis list
  enumerate_basis_default_args <- list(
    max_degree = ifelse(ncol(X) >= 20, 2, 3),
    smoothness_orders = rep(1, ncol(X)),
    # num_knots = 25,
    num_knots = 20
  )
  enumerate_basis_args <- modifyList(
    enumerate_basis_default_args,
    enumerate_basis_args
  )
  enumerate_basis_args$x <- X
  basis_list <- do.call(enumerate_basis, enumerate_basis_args)

  # basis_list <- enumerate_basis(x = X,
  #                               max_degree = ifelse(ncol(X) >= 20, 2, 3),
  #                               smoothness_orders = 1,
  #                               num_knots = hal9001:::num_knots_generator(
  #                                 max_degree = ifelse(ncol(X) >= 20, 2, 3),
  #                                 smoothness_orders = 1,
  #                                 base_num_knots_0 = 200,
  #                                 base_num_knots_1 = 20))
  #                                 #base_num_knots_1 = 50))
  # penalty_factor <- rep(1, length(basis_list))

  if (!is.null(X_weak_penalized)) {
    # NOT USED RN
    # make basis list for weakly penalized terms
    # append to basis list
    X_min <- apply(X_weak_penalized, 2, min)
    basis_list_main_terms <- enumerate_basis(
      x = X_min,
      max_degree = ifelse(
        ncol(X) >= 20, 2, 3
      ),
      smoothness_orders = 1
    )
    basis_list <- c(basis_list, basis_list_main_terms)
    penalty_factor <- c(penalty_factor, rep(
      X_weak_penalized_level,
      length(basis_list_main_terms)
    ))
  }

  # fit HAL
  fit_hal_args$X <- X
  fit_hal_args$Y <- Y
  fit_hal_args$weights <- weights
  fit_hal_args$family <- family
  fit_hal_args$basis_list <- basis_list
  fit_hal_args$return_x_basis <- TRUE
  fit_hal_args$X_unpenalized <- X_unpenalized
  hal_fit <- do.call(fit_hal, fit_hal_args)

  if (!is.null(X_unpenalized)) {
    # drop unpenalized terms from basis list, adjust indices accordingly
    hal_bases_idx <- 2:(length(hal_fit$coefs) - ncol(X_unpenalized))
    hal_basis_list <- hal_fit$basis_list[hal_fit$coefs[hal_bases_idx] != 0]
  } else {
    hal_basis_list <- hal_fit$basis_list[hal_fit$coefs[-1] != 0]
  }

  selected_idx <- which(hal_fit$coefs != 0)
  selected_unpenalized_idx <- which(
    hal_fit$coefs[length(hal_fit$basis_list) + 1:length(hal_fit$coefs)] != 0
  ) - 1
  x_basis <- cbind(1, hal_fit$x_basis)[, selected_idx, drop = FALSE]
  beta <- NULL

  if (relaxed) {
    # relaxed HAL fit
    if (family == "binomial") {
      hal_relaxed_fit <- glm.fit(x = x_basis, y = Y, family = binomial(), weights = weights, intercept = FALSE)
      beta <- coef(hal_relaxed_fit)

      # drop NA columns, collinearity
      good_idx <- which(!is.na(beta))
      beta <- beta[good_idx]
      x_basis <- x_basis[, good_idx, drop = FALSE]
      hal_basis_list <- hal_basis_list[good_idx[-1] - 1]
      pred <- x_basis %*% beta
      pred <- as.vector(1 / (1 + exp(-pred)))
    } else if (family == "gaussian") {
      hal_relaxed_fit <- glm.fit(
        x = x_basis, y = Y, family = gaussian(),
        weights = weights, intercept = FALSE
      )
      beta <- coef(hal_relaxed_fit)

      # drop NA columns, collinearity
      selected_hal_idx_good <- as.numeric(which(
        !is.na(beta[-1][seq(beta) <= length(hal_basis_list)])
      ))
      selected_unpenalized_idx_good <- as.numeric(which(
        !is.na(beta[-1][seq(beta) > length(hal_basis_list)])
      ))
      good_idx <- as.numeric(which(!is.na(beta)))
      beta <- beta[good_idx]
      x_basis <- x_basis[, good_idx, drop = FALSE]
      hal_basis_list <- hal_basis_list[selected_hal_idx_good]
      pred <- as.vector(x_basis %*% beta)
    }
  } else {
    # regular HAL fit
    beta <- hal_fit$coefs[selected_idx]
    pred <- predict(hal_fit,
      new_data = X,
      new_X_unpenalized = X_unpenalized,
      type = "response"
    )
  }

  print("number of non-zero coefficients: " %+% length(beta))

  return(list(
    beta = beta,
    hal_basis_list = hal_basis_list,
    pred = pred,
    x_basis = as.matrix(x_basis),
    selected_unpenalized_idx = selected_unpenalized_idx,
    hal_fit = hal_fit
  ))
}
