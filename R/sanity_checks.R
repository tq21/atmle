#' @title Sanity checks for input data and nodes
#'
#' @keywords internal
#'
#' @param data A \code{data.frame} containing baseline covariates \eqn{W}, binary
#' treatment indicator \eqn{A} (\eqn{A=1} for active treatment),
#' outcome \eqn{Y}, and binary indicator of whether the observation is from the
#' randomized controlled trial \eqn{S=1} or from the external data \eqn{S=0}.
#' @param S_node The column indices of the \eqn{S} node in \code{data}.
#' @param W_node The column indices of the \eqn{W} node in \code{data}.
#' @param A_node The column indices of the \eqn{A} node in \code{data}.
#' @param Y_node The column indices of the \eqn{Y} node in \code{data}.
check_data_and_nodes <- function(data, S_node, W_node, A_node, Y_node) {
  # nodes must be numeric
  if (!is.numeric(S_node)) stop("Study indicator index (S_node) must be numeric.")
  if (!is.numeric(W_node)) stop("Covariate node index (W_node) must be numeric.")
  if (!is.numeric(A_node)) stop("Treatment node index (A_node) must be numeric.")
  if (!is.numeric(Y_node)) stop("Outcome node index (Y_node) must be numeric.")

  # nodes length valid
  if (length(S_node) != 1) stop("Study indicator index (S_node) must be a single number.")
  if (length(A_node) != 1) stop("Treatment node index (A_node) must be a single number.")
  if (length(Y_node) != 1) stop("Outcome node index (Y_node) must be a single number.")

  # nodes exist in data
  if (S_node < 1 | S_node > ncol(data)) stop("Study indicator index (S_node) must be between 1 and the number of columns in data.")
  if (A_node < 1 | A_node > ncol(data)) stop("Treatment node index (A_node) must be between 1 and the number of columns in data.")
  if (Y_node < 1 | Y_node > ncol(data)) stop("Outcome node index (Y_node) must be between 1 and the number of columns in data.")
  if (min(W_node) < 1 | max(W_node) > ncol(data)) stop("Covariate node indices (W_node) must be between 1 and the number of columns in data.")

  # nodes are not duplicated
  if (length(c(S_node, W_node, A_node, Y_node)) != length(unique(c(S_node, W_node, A_node, Y_node)))) stop("Nodes must be unique.")

  # data must be a data frame
  if (!is.data.frame(data)) stop("Data must be a data frame.")
}

#' @title Sanity checks for learners
#'
#' @keywords internal
#'
#' @param theta_method The method to estimate the nuisance function
#' \eqn{\theta(W,A)=\mathbb{E}(Y\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' Default is \code{"glmnet"}.
#' @param Pi_method The method to estimate the nuisance function
#' \eqn{\Pi(S\mid W,A)=\mathbb{P}(S=1\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' @param g_method The method to estimate the nuisance function
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' @param theta_tilde_method The method to estimate the nuisance function
#' \eqn{\tilde{\theta}(W,A)=\mathbb{E}(Y\mid W,A,S=1)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' @param Q_method The method to estimate the nuisance function
#' \eqn{Q(A,W)=\mathbb{E}(Y\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a list of \code{sl3} learners for super learner-based estimation.
#' @param bias_working_model The working model for the bias estimand.
#' Either \code{"glmnet"} for lasso-based working model or \code{"HAL"} for highly adaptive
#' lasso-based working model. Default is \code{"glmnet"}.
#' @param pooled_working_model The working model for the pooled-ATE estimand.
#' Either \code{"glmnet"} for lasso-based working model or \code{"HAL"} for highly adaptive
#' lasso-based working model. Default is \code{"glmnet"}.
check_learners <- function(theta_method,
                           Pi_method,
                           g_method,
                           theta_tilde_method,
                           Q_method,
                           bias_working_model,
                           pooled_working_model) {
  # each learning method must be "glm," "glmnet," or a list of sl3 learners
  if (!is.character(theta_method) & !is.list(theta_method)) stop("theta_method must be 'glm', 'glmnet', or a list of sl3 learners.")
  if (!is.character(Pi_method) & !is.list(Pi_method)) stop("Pi_method must be 'glm', 'glmnet', or a list of sl3 learners.")
  if (!is.character(g_method) & !is.list(g_method)) stop("g_method must be 'glm', 'glmnet', or a list of sl3 learners.")
  if (!is.character(theta_tilde_method) & !is.list(theta_tilde_method)) stop("theta_tilde_method must be 'glm', 'glmnet', or a list of sl3 learners.")
  if (!is.character(Q_method) & !is.list(Q_method)) stop("Q_method must be 'glm', 'glmnet', or a list of sl3 learners.")

  # working models must be either "glmnet" or "HAL"
  if (!is.character(bias_working_model)) stop("bias_working_model must be a string.")
  if (bias_working_model != "glmnet" & bias_working_model != "HAL") stop("bias_working_model must be either 'glmnet' or 'HAL'.")
  if (!is.character(pooled_working_model)) stop("pooled_working_model must be a string.")
  if (pooled_working_model != "glmnet" & pooled_working_model != "HAL") stop("pooled_working_model must be either 'glmnet' or 'HAL'.")
}

#' Sanity checks for other arguments
check_args <- function(controls_only,
                       atmle_pooled,
                       var_method,
                       g_rct,
                       verbose) {
  # controls_only must be logical
  if (!is.logical(controls_only)) stop("controls_only must be logical.")

  # atmle_pooled must be logical
  if (!is.logical(atmle_pooled)) stop("atmle_pooled must be logical.")

  # var_method must be a string, either "ic" or "bootstrap"
  if (!is.character(var_method)) stop("var_method must be a string.")
  if (var_method != "ic" & var_method != "bootstrap") stop("var_method must be either 'ic' or 'bootstrap'.")

  # g_rct must be a numeric between 0 and 1
  if (!is.numeric(g_rct)) stop("g_rct must be numeric.")
  if (g_rct < 0 | g_rct > 1) stop("g_rct must be between 0 and 1.")

  # verbose must be logical
  if (!is.logical(verbose)) stop("verbose must be logical.")
}
