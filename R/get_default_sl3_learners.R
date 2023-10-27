#' @title Get default \code{sl3} learners
#'
#' @export
#'
#' @importFrom sl3 Lrnr_earth Lrnr_xgboost
#'
#' @param family A character string specifying the family of the outcome.
get_default_sl3_learners <- function(family) {
  learner_list <- list(
    #Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0),
    Lrnr_earth$new(degree = 3),
    Lrnr_gam$new(),
    Lrnr_glm$new())

  return(learner_list)
}
