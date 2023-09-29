options(sl3.verbose = TRUE)
source("utils.R")
set.seed(29857)

data <- generate_data(1000, 1.5, 2.6, 0.67)
#data <- generate_realistic_data(1.5, n_rct = 500, n_rwd = 2000, g_rct = 0.67, "param_simple")
S_node = 1
W_node = c(2, 3, 4, 5)
A_node = 6
Y_node = 7
nuisance_method="glm"
working_model="lasso"
p_rct=0.67
verbose=TRUE
transform=TRUE

#source("utils_positivity.R")

B <- 100
n <- 2000
bias <- 0.7
nuisance_method = "sl3"
working_model = "lasso"
pRCT = 0.67
verbose = TRUE

tmp <- run_sim(B = B,
               n = n,
               bA = 1.5,
               bias = bias,
               nuisance_method = nuisance_method,
               working_model = working_model,
               pRCT = pRCT,
               verbose = verbose,
               method = "atmle")
mean(tmp$psi_coverage)
var(tmp$psi_est)
hist(tmp$psi_est)
mean(tmp$psi_est)-1.5
var(tmp$psi_est)+(mean(tmp$psi_est)-1.5)^2

tmp_2 <- run_sim(B = B,
                 n = n,
                 bA = 1.5,
                 bias = bias,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 pRCT = pRCT,
                 verbose = verbose,
                 method = "escvtmle")
mean(tmp_2$escvtmle_prop_selected)
mean(tmp_2$psi_coverage)
var(tmp_2$psi_est)
hist(tmp_2$psi_est)
mean(tmp_2$psi_est)-1.5
var(tmp_2$psi_est)+(mean(tmp_2$psi_est)-1.5)^2

tmp_3 <- run_sim(B = B,
                 n = n,
                 bA = 1.5,
                 bias = bias,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 pRCT = pRCT,
                 verbose = verbose,
                 method = "atmle_tmle")
mean(tmp_3$psi_coverage)
var(tmp_3$psi_est)
hist(tmp_3$psi_est)
mean(tmp_3$psi_est)-1.5
var(tmp_3$psi_est)+(mean(tmp_3$psi_est)-1.5)^2

B <- 200
covered <- vector(length = B)
for (i in 1:B) {
  data <- generate_realistic_data(1.5, n_rct = 500, n_rwd = 2000, g_rct = 0.67, 0)
  S_node = 1
  W_node = c(2, 3, 4, 5)
  A_node = 6
  Y_node = 7
  p_rct=0.67
  res <- atmle(data = data,
               S_node = S_node,
               W_node = W_node,
               A_node = A_node,
               Y_node = Y_node,
               nuisance_method="glm",
               working_model="lasso",
               p_rct=0.67,
               verbose=FALSE,
               transform=TRUE)
  if (res$lower <= 1.5 & res$upper >= 1.5) {
    print("covered")
    covered[i] <- 1
  } else {
    print("not covered")
    covered[i] <- 0
  }
}
mean(covered)

