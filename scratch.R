source("utils.R")
set.seed(29857)

data <- generate_data(2000, 1.5, 0, 0.5)
S_node = 1
W_node = c(2, 3, 4, 5)
A_node = 6
Y_node = 7
nuisance_method="glm"
working_model="lasso"
p_rct=0.5
verbose=TRUE
transform=TRUE

source("utils_positivity.R")

tmp <- run_sim(B = 100,
               n = 1000,
               bA = 1.5,
               bias = 0,
               nuisance_method = "glm",
               working_model = "lasso",
               pRCT = 0.5,
               verbose = TRUE,
               method = "atmle")
mean(tmp$psi_coverage)
var(tmp$psi_est)
hist(tmp$psi_est)
mean(tmp$psi_est)-1.5
var(tmp$psi_est)+(mean(tmp$psi_est)-1.5)^2

tmp_2 <- run_sim(B = 100,
                 n = 1000,
                 bA = 1.5,
                 bias = 0,
                 nuisance_method = "glm",
                 working_model = "lasso",
                 pRCT = 0.5,
                 verbose = TRUE,
                 method = "escvtmle")
mean(tmp_2$escvtmle_prop_selected)
mean(tmp_2$psi_coverage)
var(tmp_2$psi_est)
hist(tmp_2$psi_est)
mean(tmp_2$psi_est)-1.5
var(tmp_2$psi_est)+(mean(tmp_2$psi_est)-1.5)^2

tmp_3 <- run_sim(B = 100,
                 n = 500,
                 bA = 1.5,
                 bias = 0,
                 nuisance_method = "glm",
                 working_model = "lasso",
                 pRCT = 0.67,
                 verbose = TRUE,
                 method = "atmle_tmle")
mean(tmp_3$psi_coverage)
var(tmp_3$psi_est)
hist(tmp_3$psi_est)
mean(tmp_3$psi_est)-1.5
var(tmp_3$psi_est)+(mean(tmp_3$psi_est)-1.5)^2

tmp_4 <- run_sim(B = 100,
                 n = 500,
                 bA = 1.5,
                 bias = 1.8,
                 pRCT = 0.5,
                 nuisance_method = "glm",
                 working_model = "lasso",
                 verbose = TRUE,
                 method = "rct_only")
mean(tmp_4$psi_coverage)
var(tmp_4$psi_est)
hist(tmp_4$psi_est)
mean(tmp_4$psi_est)-1.5
var(tmp_4$psi_est)+(mean(tmp_4$psi_est)-1.5)^2

tmp_5 <- run_sim(B = 100,
                 n = 1000,
                 bA = 1.5,
                 bias = 1.8,
                 nuisance_method = "glm",
                 working_model = "lasso",
                 pRCT = 0.5,
                 verbose = TRUE,
                 method = "tmle")
mean(tmp_5$psi_coverage)
var(tmp_5$psi_est)
hist(tmp_5$psi_est)
mean(tmp_5$psi_est)-1.5
var(tmp_5$psi_est)+(mean(tmp_5$psi_est)-1.5)^2



pound_est <- vector(length = 100)
tilde_est <- vector(length = 100)
for (i in 1:100) {
  print(i)
  data <- generate_data(500, 1.5, 10, 0.67)
  res <- atmle(data = data,
               S_node = 1,
               W_node = c(2, 3, 4, 5),
               A_node = 6,
               Y_node = 7,
               nuisance_method="glm",
               working_model="lasso",
               p_rct=0.67,
               verbose=FALSE)
  pound_est[i] <- res$psi_pound_est
  tilde_est[i] <- res$psi_tilde_est
}
hist(pound_est)
hist(tilde_est)

data_n <- data.frame(n = seq(200, 5000, 50),
                     tmp_bias = unlist(tmp_converge$all_psi_est))
truth <- 1.5

ggplot(data_n, aes(x = n, y = tmp_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = truth, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "est") +
  theme_minimal() +
  theme(text = element_text(size = 16))




tmp <- run_sim(B = 1,
               n = 20000,
               bA = 1.5,
               bias = "parametric",
               nuisance_method = "glm",
               working_model = "glm",
               verbose=TRUE)

tmp <- run_sim(B = 200,
               n = 500,
               bA = 1.5,
               bias = "complex",
               nuisance_method = "glm",
               working_model = "lasso",
               verbose=TRUE)



data <- generate_data(500, 1.5, "complex")
S_node = 1
W_node = c(2, 3, 4, 5)
A_node = 6
Y_node = 7
nuisance_method="lasso"
working_model="lasso"
p_rct=0.5
verbose=TRUE

rct_only(data,
         S_node,
         W_node,
         A_node,
         Y_node,
         nuisance_method="glm",
         working_model="lasso",
         p_rct=0.5,
         verbose=TRUE)

atmle(data,
      S_node,
      W_node,
      A_node,
      Y_node,
      nuisance_method="glm",
      working_model="lasso",
      p_rct=0.5,
      verbose=TRUE)


library(rlearner)
n = 100
p = 10

rlasso_fit = rlasso(as.matrix(W), A, Y, p_hat = g, m_hat = theta_tilde)
rlasso_est = predict(rlasso_fit, as.matrix(W))
mean(rlasso_est)

rlasso_fit = rlasso(as.matrix(W), A, Y)
rlasso_est = predict(rlasso_fit, as.matrix(W))
mean(rlasso_est)


rboost_fit = rboost(as.matrix(W), A, Y)
rboost_est = predict(rboost_fit, as.matrix(W))
mean(rboost_est)




library(origami)
library(purrr)

# cross-fit
folds <- make_folds(n = n, V = 5, strata_ids = A)
theta_tilde_all <- map_dfr(folds, function(.x) {
  train_idx <- .x$training_set
  valid_idx <- .x$validation_set

  # use super learner
  Q_lib <- list(
    mean = make_learner(Lrnr_mean),
    xgb = make_learner(Lrnr_xgboost, nrounds = 20, maxdepth = 6),
    ranger_small = make_learner(Lrnr_ranger, num.trees = 500),
    lasso_fast = make_learner(Lrnr_glmnet, nfold = 3)
  )
  Q_lrnr <- Lrnr_sl$new(learners = Q_lib)
  task <- make_sl3_Task(data[train_idx,], c("W1", "W2", "W3", "W4", "A"), "Y")
  task_valid <- make_sl3_Task(data[valid_idx,], c("W1", "W2", "W3", "W4", "A"), "Y")
  sl_fit <- Q_lrnr$train(task)
  theta_tilde_sl <- sl_fit$predict(task_valid)

  #theta_tilde_fit <- fit_hal(X = W[train_idx,], Y = Y[train_idx], family = "gaussian", smoothness_orders = 0)
  #theta_tilde_pred <- predict(theta_tilde_fit, new_data = W[valid_idx,])

  return(as.data.frame(cbind(valid_idx, theta_tilde_sl)))
})
theta_tilde_all <- theta_tilde_all[order(theta_tilde_all$valid_idx),]

psi <- learn_psi_tilde(W, A, Y, g, theta_tilde_all$theta_tilde_pred)
psi_est <- mean(psi$pred)


psi_eic <- get_eic_psi_tilde(psi, g, theta_tilde, Y, A, n)
psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
psi_ci_lower <- psi_est-1.96*psi_se
psi_ci_upper <- psi_est+1.96*psi_se




# use super learner
Q_lib <- list(
  mean = make_learner(Lrnr_mean),
  glm = make_learner(Lrnr_glm_fast),
  xgb = make_learner(Lrnr_xgboost, nrounds = 20, maxdepth = 6),
  ranger_small = make_learner(Lrnr_ranger, num.trees = 500),
  lasso_fast = make_learner(Lrnr_glmnet, nfold = 3)
)
Q_lrnr <- Lrnr_sl$new(learners = Q_lib)
task <- make_sl3_Task(data, c("W1", "W2", "W3", "W4", "A"), "Y")
sl_fit <- Q_lrnr$train(task)
theta_tilde_sl <- sl_fit$predict()
sqrt(mean((Y-theta_tilde_sl)^2))

psi_tilde_sl <- learn_psi_tilde(W, A, Y, g, theta_tilde_sl)
mean(psi_tilde$pred)
mean(psi_tilde_sl$pred)



no_bias_n_increase <- run_sim_n_increase(B = 1,
                                         n_min = 200,
                                         n_max = 5000,
                                         n_step = 50,
                                         bA = 1.5,
                                         bias = 0,
                                         nuisance_method = "lasso",
                                         working_model = "lasso",
                                         verbose=TRUE)

constant_bias_n_increase <- run_sim_n_increase(B = 1,
                                               n_min = 200,
                                               n_max = 5000,
                                               n_step = 50,
                                               bA = 1.5,
                                               bias = 1.8,
                                               nuisance_method = "lasso",
                                               working_model = "lasso",
                                               verbose=TRUE)

data_n <- data.frame(n = seq(200, 5000, 50),
                     no_bias = unlist(no_bias_n_increase$all_psi_est),
                     constant_bias = unlist(constant_bias_n_increase$all_psi_est))

p_no_bias_bias_converge <- ggplot(data_n, aes(x = n, y = no_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = truth, color = "red", linetype = "dashed", linewidth = 1) +
  #scale_y_continuous(breaks = seq(1.2, 1.8, 0.1), limits = c(1.2, 1.8)) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_constant_bias_bias_converge <- ggplot(data_n, aes(x = n, y = constant_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = truth, color = "red", linetype = "dashed", linewidth = 1) +
  #scale_y_continuous(breaks = seq(1.2, 1.8, 0.1), limits = c(1.2, 1.8)) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))


set.seed(123)
W <- runif(1000)
A <- rbinom(1000, 1, 0.6*W)
Y <- 0.5+1.3*A+0.7*W
fit <- lm(A ~ W)
coef(fit)


