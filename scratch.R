options(sl3.verbose = TRUE)
source("utils.R")
set.seed(29857)

data <- generate_realistic_data(1.5, n_rct = 200, n_rwd = 2000, g_rct = 0.67, bias = "param_simple", TRUE)
S_node = 1
W_node = c(2, 3, 4, 5)
A_node = 6
Y_node = 7
nuisance_method="glm"
working_model="glmnet"
g_rct=0.67
verbose=TRUE
transform=TRUE
controls_only = TRUE

#source("utils_positivity.R")

B <- 200
n_rct <- 200
n_rwd <- 1000
ate <- 1.5
bias <- "param_complex"
nuisance_method = "glm"
working_model = "glmnet"
g_rct = 0.67
verbose = TRUE
controls_only = TRUE

tmp <- run_sim(B = B,
               n_rct = n_rct,
               n_rwd = n_rwd,
               ate = ate,
               bias = bias,
               nuisance_method = nuisance_method,
               working_model = working_model,
               g_rct = g_rct,
               controls_only = controls_only,
               verbose = verbose,
               method = "atmle")
mean(tmp$psi_coverage)
var(tmp$psi_est)
hist(tmp$psi_est)
mean(tmp$psi_est)-1.5
var(tmp$psi_est)+(mean(tmp$psi_est)-1.5)^2

tmp_2 <- run_sim(B = B,
                 n_rct = n_rct,
                 n_rwd = n_rwd,
                 ate = ate,
                 bias = bias,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 g_rct = g_rct,
                 controls_only = controls_only,
                 verbose = verbose,
                 method = "escvtmle")
mean(tmp_2$escvtmle_prop_selected)
mean(tmp_2$psi_coverage)
var(tmp_2$psi_est)
hist(tmp_2$psi_est)
mean(tmp_2$psi_est)-1.5
var(tmp_2$psi_est)+(mean(tmp_2$psi_est)-1.5)^2

tmp_3 <- run_sim(B = B,
                 n_rct = n_rct,
                 n_rwd = n_rwd,
                 ate = ate,
                 bias = bias,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 g_rct = g_rct,
                 verbose = verbose,
                 method = "nonparametric")
mean(tmp_3$psi_coverage)
var(tmp_3$psi_est)
hist(tmp_3$psi_est)
mean(tmp_3$psi_est)-1.5
var(tmp_3$psi_est)+(mean(tmp_3$psi_est)-1.5)^2

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

# real data comparison
data(wash)
#For unbiased external controls, use:
dat <- wash[which(wash$study %in% c(1,2)),]
dat$study[which(dat$study==2)]<-0
dat$study[which(dat$study==3)]<-0
dat$sex <- as.numeric(dat$sex)-1
dat$momedu <- as.numeric(dat$momedu)-1
dat$hfiacat <- as.numeric(dat$hfiacat)-1

set.seed(2022)
res_escvtmle <- ES.cvtmle(txinrwd=TRUE,
                          data=dat, study="study",
                          covariates=c("aged", "sex", "momedu", "hfiacat"),
                          treatment_var="intervention", treatment=1,
                          outcome="laz", NCO="Nlt18scale",
                          Delta=NULL, Delta_NCO=NULL,
                          pRCT=0.5, V=5, Q.SL.library=c("SL.glm"),
                          g.SL.library=c("SL.glm"), Q.discreteSL=TRUE, g.discreteSL=TRUE,
                          family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                          comparisons = list(c(1),c(1,0)), adjustnco = FALSE, target.gwt = TRUE)
print.EScvtmle(res_escvtmle)

res_atmle <- atmle(data = dat,
                   S_node = 2,
                   W_node = c(4, 5, 6, 7),
                   A_node = 1,
                   Y_node = 3,
                   nuisance_method="glm",
                   working_model="lasso",
                   g_rct=0.5,
                   verbose=TRUE)
print("ATMLE estimate: " %+% round(res_atmle$est, 3) %+% " (" %+%
        round(res_atmle$lower, 3) %+% ", " %+% round(res_atmle$upper, 3) %+% ")")

res_atmle$est
res_atmle$lower
res_atmle$upper
