library(devtools)
library(EScvtmle)
library(sl3)
load_all()
#set.seed(123)
sim_data <- function(ate,
                     n,
                     rct,
                     g_rct,
                     bias,
                     controls_only) {
  # error
  UY <- rnorm(n, 0, 1)

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)

  # study indicator S and treatment A
  if (rct) {
    S <- rep(1, n)
    A <- rbinom(n, 1, g_rct)
  } else {
    S <- rep(0, n)
    if (controls_only) {
      A <- rep(0, n)
    } else {
      A <- rbinom(n, 1, plogis(0.5*W1))
    }
  }

  # bias term for RWD data
  if (bias == "a") {
    b <- 0.2+0.2*W1*(1-A)
  } else if (bias == "b") {
    b <- 0.5+3.1*W1*(1-A)+0.8*W3
  }

  obs_prob <- plogis(2)#plogis(2+0.5*W1)
  delta <- rbinom(n, 1, obs_prob)

  # outcome
  Y <- 2.5+0.9*W1+1.1*W2+2.7*W3+ate*A+UY+(1-S)*b
  #Y[delta == 0] <- NA

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

controls_only <- TRUE

data_rct <- sim_data(ate = 1.5,
                     n = 1000,
                     rct = TRUE,
                     g_rct = 0.67,
                     bias = "a",
                     controls_only = controls_only)
data_rwd <- sim_data(ate = 1.5,
                     n = 5000,
                     rct = FALSE,
                     g_rct = 0.67,
                     bias = "a",
                     controls_only = controls_only)
data <- rbind(data_rct, data_rwd)

res <- atmle_new(data = data,
                 S = "S",
                 W = c("W1", "W2", "W3"),
                 A = "A",
                 Y = "Y",
                 controls_only = controls_only,
                 family = "gaussian",
                 bias_working_model = "glmnet",
                 pooled_working_model = "glmnet",
                 target_gwt = TRUE,
                 verbose = FALSE,
                 target_method = "oneshot",
                 browse = FALSE)

res_escvtmle <- ES.cvtmle(txinrwd = !controls_only,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2", "W3"),
                          treatment_var = "A",
                          treatment = 1,
                          outcome = "Y",
                          pRCT = 0.67,
                          family = "gaussian",
                          Q.SL.library = c("SL.glm"),
                          g.SL.library = c("SL.glm"),
                          d.SL.library.RCT = c("SL.glm"),
                          d.SL.library.RWD = c("SL.glm"),
                          Q.discreteSL = TRUE,
                          g.discreteSL = TRUE,
                          d.discreteSL = TRUE,
                          V = 5)

res$upper-res$lower
as.numeric(res_escvtmle$CI$b2v[2]-res_escvtmle$CI$b2v[1])
res$est
res_escvtmle

library(devtools)
load_all()
sim_data <- function(n,
                     counter_A = NULL) {
  # error
  UY <- rnorm(n, 0, 0.2)

  # baseline covariates
  W1 <- round(runif(n, -1, 1), 3)
  W2 <- round(runif(n, -1, 1), 3)
  W3 <- round(runif(n, -1, 1), 3)

  # treatment
  if (is.null(counter_A)) {
    A <- rbinom(n, 1, plogis(-0.25*W1+0.5*W2))
  } else {
    A <- rep(counter_A, n)
  }

  # outcome
  cate <- 0.5*as.numeric(W1 >= 0.2)*(W1-0.2)+
    1.1*as.numeric(W2 >= 0.5)*(W2-0.5)+
    2.1*as.numeric(W1 >= -0.5)*(W2+0.5)*as.numeric(W3 >= 0)*W3+
    1.4*as.numeric(W2 >= 0.9)*(W2-0.9)*as.numeric(W3 >= 0.1)*(W3-0.1)
  Y <- 1.9-0.4*A+0.9*W1+1.4*W2+2.1*W3+A*cate+UY

  data <- data.frame(W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

get_truth <- function() {
  data_A1 <- sim_data(1e7, counter_A = 1)
  data_A0 <- sim_data(1e7, counter_A = 0)
  return(mean(data_A1$Y- data_A0$Y))
}

theta_lasso <- function(W,
                        Y,
                        foldid,
                        family) {

  # fit theta=E(Y|W) and obtain cross-fitted predictions
  theta_fit <- cv.glmnet(x = as.matrix(W),
                         y = Y,
                         foldid = foldid,
                         keep = TRUE,
                         alpha = 1,
                         family = family)
  non_na_idx <- !is.na(colSums(theta_fit$fit.preval))
  lambda_min <- theta_fit$lambda[which.min(theta_fit$cvm[non_na_idx])]
  pred <- theta_fit$fit.preval[, non_na_idx][, theta_fit$lambda[non_na_idx] == lambda_min]
  if (family == "binomial") {
    pred <- plogis(pred)
  }

  return(pred)
}

g_lasso <- function(W,
                    A,
                    foldid) {

  # fit g1W=P(A=1|W) and obtain cross-fitted predictions
  g_fit <- cv.glmnet(x = as.matrix(W),
                     y = A,
                     foldid = foldid,
                     family = "binomial",
                     keep = TRUE,
                     alpha = 1)
  non_na_idx <- !is.na(colSums(g_fit$fit.preval))
  lambda_min <- g_fit$lambda[which.min(g_fit$cvm[non_na_idx])]
  pred <- g_fit$fit.preval[, non_na_idx][, g_fit$lambda[non_na_idx] == lambda_min]
  pred <- plogis(pred)

  return(pred)
}


truth <- get_truth()

library(origami)
library(purrr)
library(glmnet)
library(hal9001)
library(tmle)
set.seed(123)
n <- 5000
data <- sim_data(n)
W <- data[, grep("W", colnames(data))]
A <- data$A
Y <- data$Y




# make folds
folds <- make_folds(n = n, V = 10, strata_ids = A)
foldid <- unlist(map(folds, function(.fold) {
  rep(.fold$v, length(.fold$validation_set))
}))
idx <- unlist(map(folds, function(.fold) {
  .fold$validation_set
}))
foldid <- foldid[idx]

theta <- theta_lasso(W = W,
                     Y = Y,
                     foldid = foldid,
                     family = "gaussian")
g1W <- g_lasso(W = W,
               A = A,
               foldid = foldid)
res <- learn_tau_A(W = W,
                   A = A,
                   Y = Y,
                   g1W = g1W,
                   delta = rep(1, n),
                   theta = theta,
                   method = "HAL",
                   foldid = foldid,
                   weights = rep(1, n),
                   enumerate_basis_args = list(max_degree = 2,
                                               smoothness_orders = 1,
                                               num_knots = c(20, 5)),
                   pooled_working_model_formula = NULL,
                   target_method = "oneshot",
                   eic_method = "svd_pseudo_inv",
                   verbose = FALSE,
                   browse = FALSE)
res_relax <- learn_tau_A(W = W,
                         A = A,
                         Y = Y,
                         g1W = g1W,
                         delta = rep(1, n),
                         theta = theta,
                         method = "HAL",
                         foldid = foldid,
                         weights = rep(1, n),
                         enumerate_basis_args = list(max_degree = 2,
                                                     smoothness_orders = 1,
                                                     num_knots = c(20, 5)),
                         pooled_working_model_formula = NULL,
                         target_method = "relaxed",
                         eic_method = "svd_pseudo_inv",
                         verbose = FALSE,
                         browse = FALSE)
res_tmle <- tmle(W = W, A = A, Y = Y, family = "gaussian",
                 Q.SL.library = c("SL.glm", "SL.xgboost"), g.SL.library = c("SL.glm", "SL.xgboost"))
res_tmle$estimates$ATE$psi

a <- proc.time()
tau <- rHAL(W = W,
            A = A,
            Y = Y,
            g1W = g1W,
            theta = theta,
            foldid = foldid,
            use_weight = TRUE)
b <- proc.time()
tau_weight <- rHAL(W = W,
                   A = A,
                   Y = Y,
                   g1W = g1W,
                   theta = theta,
                   foldid = foldid,
                   use_weight = FALSE)
c <- proc.time()
mean(tau$pred)
mean(tau_weight$pred)

predict((tau), newx = W[1:10,])
predict((tau_weight), newx = W[1:10,])
