#' @title Get default \code{sl3} learners
#'
#' @description Function to get a default super learner library for a given
#' family of the outcome.
#'
#' @keywords internal
#'
#' @importFrom sl3 Lrnr_glm Lrnr_dbarts Lrnr_xgboost
#'
#' @param family A character string specifying the family of the outcome.
#' Either \code{"gaussian"} or \code{"binomial"}.
get_default_sl3_learners <- function(family) {
  learner_list <- list(
    Lrnr_glm$new(),
    Lrnr_dbarts$new(),
    Lrnr_xgboost$new()
  )

  return(learner_list)
}
