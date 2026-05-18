#' @noRd
fit_regression <- function(data,
                           method,
                           folds,
                           covariate_nodes,
                           outcome_node,
                           weight_node = NULL,
                           subset = seq(nrow(data)),
                           bound = NULL) {

  if (!is.list(method) || is.null(method$learners)) {
    stop("`method` must be a list with a `learners` element.")
  }

  if (length(method$learners) == 1) {
    lrnr <- method$learners[[1]]
    if (!is.null(bound)) {
      lrnr_bound <- sl3::Lrnr_bound$new(bound)
      lrnr <- sl3::Pipeline$new(lrnr, lrnr_bound)
    }
  } else {
    lrnr_stack <- sl3::Stack$new(method$learners)
    if (!is.null(bound)) {
      lrnr_bound <- sl3::Lrnr_bound$new(bound)
      lrnr_stack <- sl3::Pipeline$new(lrnr_stack, lrnr_bound)
    }
    lrnr <- sl3::make_learner(
      sl3::Pipeline,
      sl3::Lrnr_cv$new(lrnr_stack),
      method$metalearner
    )
  }

  suppressWarnings({
    task <- sl3::sl3_Task$new(
      data = data[subset, , drop = FALSE],
      covariates = covariate_nodes,
      outcome = outcome_node,
      weights = weight_node,
      folds = folds
    )
  })
  suppressMessages(fit_obj <- lrnr$train(task))

  return(invisible(fit_obj))
}
