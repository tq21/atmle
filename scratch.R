options(sl3.verbose = TRUE)
source("utils.R")
set.seed(29857)

large_n <- 10000000
W1 <- rnorm(large_n, 0, 1)
W2 <- rnorm(large_n, 0, 1)
W3 <- rnorm(large_n, 0, 1)
W4 <- rnorm(large_n, 0, 1)
ate <- mean(plogis(1.5+0.8*W1-1.1*W2+0.9*W3-1.3*W4)-plogis(-2+0.8*W1-1.1*W2+0.9*W3-1.3*W4))

data <- sim_binary_outcome(1.5, 5000, 0.2, 0.67, "param_simple", FALSE)

#data <- generate_four_covs(1.5, 50000, 0.5, g_rct = 0.67, bias = 0, FALSE, "gaussian")
S_node = 1
W_node = c(2, 3, 4, 5)
A_node = 6
Y_node = 7
nuisance_method="glm"
working_model="glmnet"
g_rct=0.67
verbose=TRUE
transform=TRUE
controls_only = FALSE
family= "binomial"

#source("utils_positivity.R")

B <- 200
n <- 2000
ate <- ate
bias <- "param_complex"
nuisance_method = "sl3"
working_model = "glmnet"
g_rct = 0.67
verbose = TRUE
controls_only = FALSE
num_covs <- 4

tmp_1 <- run_sim(B = B,
                 n = n,
                 ate = ate,
                 bias = bias,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 g_rct = g_rct,
                 controls_only = controls_only,
                 var_method = "ic",
                 num_covs = num_covs,
                 verbose = verbose,
                 method = "oracle_atmle")
mean(tmp_1$psi_coverage)
var(tmp_1$psi_est)+(mean(tmp_1$psi_est)-1.5)^2

tmp <- run_sim(B = B,
               n = n,
               ate = ate,
               bias = bias,
               nuisance_method = nuisance_method,
               working_model = working_model,
               g_rct = g_rct,
               controls_only = controls_only,
               var_method = "ic",
               num_covs = num_covs,
               verbose = verbose,
               method = "atmle")
mean(tmp$psi_coverage)
var(tmp$psi_est)
hist(tmp$psi_est)
mean(tmp$psi_est)-1.5
var(tmp$psi_est)+(mean(tmp$psi_est)-ate)^2
tmp$psi_ci_upper-tmp$psi_ci_lower

tmp_2 <- run_sim(B = B,
                 n = n,
                 ate = ate,
                 bias = bias,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 g_rct = g_rct,
                 controls_only = controls_only,
                 num_covs = num_covs,
                 verbose = verbose,
                 method = "escvtmle")
mean(tmp_2$escvtmle_prop_selected)
mean(tmp_2$psi_coverage)
var(tmp_2$psi_est)
hist(tmp_2$psi_est)
mean(tmp_2$psi_est)-1.5
var(tmp_2$psi_est)+(mean(tmp_2$psi_est)-ate)^2
tmp_2$psi_ci_upper-tmp_2$psi_ci_lower

tmp_3 <- run_sim(B = B,
                 n = n,
                 ate = ate,
                 bias = bias,
                 num_covs = num_covs,
                 controls_only = controls_only,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 g_rct = g_rct,
                 verbose = verbose,
                 method = "rct_only")
mean(tmp_3$psi_coverage)
var(tmp_3$psi_est)
hist(tmp_3$psi_est)
mean(tmp_3$psi_est)-1.5
var(tmp_3$psi_est)+(mean(tmp_3$psi_est)-1.5)^2
tmp_3$psi_ci_upper-tmp_3$psi_ci_lower

B <- 300
covered <- vector(length = B)
all_res <- vector(length = B)
for (i in 1:B) {
  data <- generate_realistic_data(1.5, n_rct = 200, n_rwd = 2000, g_rct = 0.67, bias = "param_simple")
  S_node = 1
  W_node = c(2, 3, 4, 5)
  A_node = 6
  Y_node = 7
  p_rct = 0.67
  res <- atmle(data = data,
               S_node = S_node,
               W_node = W_node,
               A_node = A_node,
               Y_node = Y_node,
               atmle_pooled = TRUE,
               var_method = "bootstrap",
               nuisance_method="glm",
               working_model="lasso",
               g_rct=0.67,
               verbose=FALSE)
  all_res[i] <- res$est
  if (res$lower <= 1.5 & res$upper >= 1.5) {
    print(i %+% ": covered")
    covered[i] <- 1
  } else {
    print(i %+% ": not covered")
    covered[i] <- 0
  }
}
mean(covered)
hist(all_res)
var(all_res)

`%+%` <- function(a, b) paste0(a, b)
# real data comparison
data(wash)
#For unbiased external controls, use:
dat <- wash[which(wash$study %in% c(1,2)),]
dat$study[which(dat$study==2)]<-0
dat$study[which(dat$study==3)]<-0
dat$sex <- as.numeric(dat$sex)-1
dat$momedu <- as.numeric(dat$momedu)-1
dat$hfiacat <- as.numeric(dat$hfiacat)-1
#dat <- subset(dat, !(study == 0 & intervention == 1))

# binary outcome
dat$laz <- as.numeric(dat$laz > -2.3)

set.seed(2022)
res_escvtmle <- ES.cvtmle(txinrwd=TRUE,
                          data=dat, study="study",
                          covariates=c("aged", "sex", "momedu", "hfiacat"),
                          treatment_var="intervention", treatment=1,
                          outcome="laz", NCO="Nlt18scale",
                          Delta=NULL, Delta_NCO=NULL,
                          pRCT=0.5, V=5, Q.SL.library=c("SL.glm"),
                          g.SL.library=c("SL.glm"), Q.discreteSL=TRUE, g.discreteSL=TRUE,
                          family="binomial", family_nco="gaussian", fluctuation = "logistic",
                          comparisons = list(c(1),c(1,0)), adjustnco = FALSE, target.gwt = TRUE)
print.EScvtmle(res_escvtmle)

res_atmle <- atmle(data = dat,
                   S_node = 2,
                   W_node = c(4, 5, 6, 7),
                   A_node = 1,
                   Y_node = 3,
                   controls_only = FALSE,
                   nuisance_method="sl3",
                   working_model="glmnet",
                   g_rct=0.5,
                   verbose=TRUE)
print("ATMLE estimate: " %+% round(res_atmle$est, 3) %+% " (" %+%
        round(res_atmle$lower, 3) %+% ", " %+% round(res_atmle$upper, 3) %+% ")")

res_atmle$est
res_atmle$lower
res_atmle$upper


library(atmle)
library(EScvtmle)
library(sl3)
library(data.table)
library(glmnet)
options(sl3.verbose = TRUE)
`%+%` <- function(a, b) paste0(a, b)
set.seed(82379)

n_rct <- 200 # RCT sample size
n_rwd <- 1000 # RWD sample size
n <- n_rct + n_rwd
g_rct <- 0.67 # RCT treatment probability
true_ate <- 1.5 # true ATE

# error
UY <- rnorm(n, 0, 0.5)
U_bias <- rnorm(n, 0, 0.1)

# baseline covariates
W1 <- rnorm(n, 0, 1)
W2 <- rnorm(n, 0, 1)
W3 <- rnorm(n, 0, 1)
W4 <- rnorm(n, 0, 1)

# study indicator, S=1 for RCT, S=0 for RWD
S <- c(rep(1, n_rct), rep(0, n_rwd))

# treatments (external data has both treated and controls)
A_rct <- rbinom(n_rct, 1, g_rct)
A_rwd <- rbinom(n_rwd, 1, plogis(0.5*W1-0.9*W2))
A <- c(A_rct, A_rwd)

# outcome
Y <- 2.1+0.8*W1+2.5*W2-3.1*W3+0.9*W4+true_ate*A+UY+
  (1-S)*(A*(2.9*W1+2.3*W2+U_bias))

# data frames combining RCT and RWD
data <- data.frame(S = S,
                   W1 = W1,
                   W2 = W2,
                   W3 = W3,
                   W4 = W4,
                   A = A,
                   Y = Y)

# run A-TMLE
res <- atmle(data = data,
             S_node = 1,
             W_node = c(2, 3, 4, 5),
             A_node = 6,
             Y_node = 7,
             controls_only = FALSE,
             atmle_pooled = TRUE,
             var_method = "ic",
             nuisance_method = "sl3",
             working_model = "glmnet",
             g_rct = g_rct,
             verbose = FALSE)
print("A-TMLE ATE estimate: " %+% round(res$est, 2) %+%
        " (" %+% round(res$lower, 2) %+% ", " %+% round(res$upper, 2) %+% ")")
print("True ATE: " %+% true_ate)

# compared to ES-CVTMLE
escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2", "W3", "W4"),
                          treatment_var = "A",
                          treatment = 1,
                          outcome = "Y",
                          pRCT = g_rct,
                          family = "gaussian",
                          Q.SL.library = c("SL.glm"),
                          g.SL.library = c("SL.glm"),
                          Q.discreteSL = TRUE,
                          g.discreteSL = TRUE,
                          V = 5)
print("ES-CVTMLE ATE estimate: " %+% round(escvtmle_res$ATE$b2v, 2) %+%
        " (" %+% round(escvtmle_res$CI$b2v[1], 2) %+%
        ", " %+% round(escvtmle_res$CI$b2v[2], 2) %+% ")")
print("True ATE: " %+% true_ate)
