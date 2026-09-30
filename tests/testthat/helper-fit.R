fixture_fits <- new.env(parent = emptyenv())

make_fixture <- function(family = "gaussian", controls_only = FALSE,
                         forced = TRUE, target_method = "tmle", n_lambda = 1L,
                         target_gwt = TRUE) {
  set.seed(135)
  n <- 300
  W1 <- rbinom(n, 1, .5)
  W2 <- rnorm(n)
  S <- rbinom(n, 1, .5)
  A <- rbinom(n, 1, .5)
  if (controls_only) A[S == 0] <- 0
  mu <- .3*W1 + .2*W2 + .5*A + .2*(1-S)*A
  Y <- if (family == "gaussian") mu + rnorm(n) else rbinom(n, 1, plogis(mu-1))
  data <- data.frame(W1, W2, S, A, Y)
  fit <- atmle_ate_fusion$new(data, "S", c("W1", "W2"), "A", "Y", family, 3)
  learner <- list(learners = list(sl3::Lrnr_glm$new()))
  suppressWarnings(fit$run(learner, learner, learner, learner,
    A_cate_args = list(max_degree = 1L, smoothness_orders = 1L, num_knots = 3L),
    S_cate_args = list(max_degree = 1L, smoothness_orders = 1L, num_knots = 3L,
                       force_A = forced),
    target_method = target_method, target_gwt = target_gwt, n_lambda = n_lambda,
    max_iter = 200L, verbose = FALSE))
  fit
}

fixture_fit <- function() {
  if (is.null(fixture_fits$fit)) fixture_fits$fit <- make_fixture()
  fixture_fits$fit$clone(deep = TRUE)
}
