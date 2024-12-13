#' @title Get default \code{sl3} learners
#'
#' @description Function to get a default super learner library for a given
#' family of the outcome.
#'
#' @keywords internal
#'
#' @importFrom sl3 Lrnr_earth Lrnr_xgboost Lrnr_gam Lrnr_glm
#'
#' @param family A character string specifying the family of the outcome.
#' Either \code{"gaussian"} or \code{"binomial"}.
#'
#' @examples
#' lrnrs <- get_default_sl3_learners("gaussian")
get_default_sl3_learners <- function(family) {
  learner_list <- list(
    Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0),
    #ranger_small = make_learner(Lrnr_ranger, num.trees = 500),
    ridge_fast = make_learner(Lrnr_glmnet, nfold = 3, alpha = 0),
    enet_fast = make_learner(Lrnr_glmnet, nfold = 3, alpha = 0.5),
    #bart = Lrnr_dbarts$new(ndpost = 1000, verbose = FALSE),
    #Lrnr_earth$new(degree = 3),
    Lrnr_gam$new()#,
    #Lrnr_glm$new()
  )

  return(learner_list)
}
