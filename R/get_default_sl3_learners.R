get_default_sl3_learners <- function(family) {
  learner_list <- list(
    Lrnr_earth$new(degree = 3, family = family),
    Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0))

  return(learner_list)
}
