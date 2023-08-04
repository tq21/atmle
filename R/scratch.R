library(hal9001)
library(sl3)
library(data.table)

set.seed(123)
S <- rbinom(100, 1, 0.5)
A <- rbinom(100, 1, 0.5)
W <- rnorm(100)
X <- matrix(c(S, A, W), nrow = 100, ncol = 3)
Y <- rnorm(100)

binomial_lrnr_stack <- Stack$new(list(Lrnr_earth$new(degree = 2, family = "binomial"),
                                 Lrnr_gam$new(family = "binomial"),
                                 Lrnr_ranger$new(),
                                 Lrnr_xgboost$new(max_depth = 4, nrounds = 20)))
gaussian_lrnr_stack <- Stack$new(list(Lrnr_earth$new(degree = 2, family = "gaussian"),
                                      Lrnr_gam$new(family = "gaussian"),
                                      Lrnr_ranger$new(),
                                      Lrnr_xgboost$new(max_depth = 4, nrounds = 20)))
sl3_lrnr_Pi <- make_learner(Pipeline,
                            Lrnr_cv$new(binomial_lrnr_stack),
                            Lrnr_cv_selector$new(loss_loglik_binomial))
sl3_lrnr_Q_bar <- make_learner(Pipeline,
                               Lrnr_cv$new(gaussian_lrnr_stack),
                               Lrnr_cv_selector$new(loss_squared_error))


# PARAMETRIC -------------------------------------------------------------------

# WORKING MODEL FOR TAU, SAME WORKING MODEL FOR PSI_TILDE ----------------------
# learn Pi(1|W,A)=P(S=1|W,A)
fit_Pi <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), S, "binomial")
Pi_pred <- fit_Pi$pred

# learn Q_bar(S,W,A)=E(Y|S,W,A)
fit_Q_bar <- fit_relaxed_hal(as.matrix(data.table(S = S, W, A = A)), Y, "gaussian")

x_basis_S1 <- make_counter_design_matrix(fit_Q_bar$basis_list, as.matrix(data.table(S = 1, W, A = A)))
x_basis_S0 <- make_counter_design_matrix(fit_Q_bar$basis_list, as.matrix(data.table(S = 0, W, A = A)))
x_basis_A1 <- make_counter_design_matrix(fit_Q_bar$basis_list, as.matrix(data.table(S = S, W, A = 1)))
x_basis_A0 <- make_counter_design_matrix(fit_Q_bar$basis_list, as.matrix(data.table(S = S, W, A = 0)))

Q_bar_S1_pred <- x_basis_S1 %*% fit_Q_bar$beta
Q_bar_S0_pred <- x_basis_S0 %*% fit_Q_bar$beta
Q_bar_A1_pred <- x_basis_A1 %*% fit_Q_bar$beta
Q_bar_A0_pred <- x_basis_A0 %*% fit_Q_bar$beta

Q_bar_S_diff_pred <- Q_bar_S1_pred - Q_bar_S0_pred
Q_bar_A_diff_pred <- Q_bar_A1_pred - Q_bar_A0_pred

# target parameter point estimate
psi_pound_pred <- mean(Pi_pred * Q_bar_S_diff_pred)
psi_tilde_pred <- mean(Q_bar_A_diff_pred)
psi_pred <- psi_tilde_pred - psi_pound_pred

# TODO: inference







# using super learner

task_Pi <- sl3_Task$new(data.table(S = S, W, A = A),
                        covariates = c(colnames(W), "A"), outcome = "S",
                        outcome_type = "binomial")
fit_Pi <- sl3_lrnr_Pi$train(task_Pi)
Pi_pred <- fit_Pi$predict(task_Pi)

# learn Q_bar(S,W,A)=E(Y|S,W,A)
task_Q_bar <- sl3_Task$new(data.table(Y = Y, S = S, W, A = A),
                           covariates = c("S", colnames(W), "A"), outcome = "Y",
                           outcome_type = "gaussian")
task_Q_bar_S1 <- sl3_Task$new(data.table(Y = Y, S = 1, W, A = A),
                              covariates = c("S", colnames(W), "A"), outcome = "Y",
                              outcome_type = "gaussian")
task_Q_bar_S0 <- sl3_Task$new(data.table(Y = Y, S = 0, W, A = A),
                              covariates = c("S", colnames(W), "A"), outcome = "Y",
                              outcome_type = "gaussian")
task_Q_bar_A1 <- sl3_Task$new(data.table(Y = Y, S = S, W, A = 1),
                              covariates = c("S", colnames(W), "A"), outcome = "Y",
                              outcome_type = "gaussian")
task_Q_bar_A0 <- sl3_Task$new(data.table(Y = Y, S = S, W, A = 0),
                              covariates = c("S", colnames(W), "A"), outcome = "Y",
                              outcome_type = "gaussian")
fit_Q_bar <- sl3_lrnr_Q_bar$train(task_Q_bar)
Q_bar_S1_pred <- fit_Q_bar$predict(task_Q_bar_S1)
Q_bar_S0_pred <- fit_Q_bar$predict(task_Q_bar_S0)
Q_bar_A1_pred <- fit_Q_bar$predict(task_Q_bar_A1)
Q_bar_A0_pred <- fit_Q_bar$predict(task_Q_bar_A0)
Q_bar_S_diff_pred <- Q_bar_S1_pred - Q_bar_S0_pred
Q_bar_A_diff_pred <- Q_bar_A1_pred - Q_bar_A0_pred

psi_pound_pred <- mean(Pi_pred * Q_bar_S_diff_pred)
psi_tilde_pred <- mean(Q_bar_A_diff_pred)
psi_pred <- psi_tilde_pred - psi_pound_pred

# evaluate EIC


