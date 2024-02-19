#devtools::install_github("rje42/causl", build_vignettes = F)
library(causl)
library(data.table)
library(atmle)
library(EScvtmle)
library(devtools)
load_all()

`%+%` <- function(a, b) paste0(a, b)

# list of coefficient for the unobserved confounder to introduce bias
coeff.U.list <- c(0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1,1.5,2,2.5)

# function to generate data
generate.data.list.HTE <- function(n.r, n.o, seed, coeff.U.list,
                                   HTE = FALSE,
                                   rct.treat.prob = 0.5,
                                   X.diff, U.diff,
                                   U1.mean.shift = 0,
                                   U2.mean.shift = 0){

  data.list <- list()

  for (i in 1:length(coeff.U.list)) {

    coeff.U <- coeff.U.list[i]

    family <- list(c(1,5), c(5,5), 1, 1)
    forms_obs_d <- forms_exp_d <- list(c(X1 ~ 1, X2 ~ 1),
                                       list(A ~ X1 + X2  + U1 + U2, X3 ~ 1),
                                       Y ~ A + X3 + A:X3 + U1 + U2,
                                       ~ 1)
    forms_exp_d[[2]][1] <- list(A ~ 1)

    pars_obs_d <- pars_exp_d <- list(X1 = list(beta = 0, phi = 1),
                                     X2 = list(beta = 0),
                                     X3 = list(beta = 0),
                                     A = list(beta = c(0,1,1,coeff.U,coeff.U)),
                                     Y = list(beta = c(0.5,0.2,1,1,1,0), phi = 1),
                                     cop = list(beta = matrix(c(0.1,2,1),nrow = 1)))
    pars_exp_d$A <- list(beta = logit(rct.treat.prob))

    if (HTE){
      pars_obs_d$Y$beta[length(pars_obs_d$Y$beta)] <- 0.1
      pars_exp_d$Y$beta[length(pars_exp_d$Y$beta)] <- 0.1
    }

    if (X.diff) {
      pars_obs_d$X1 <- list(beta = 0.5, phi = 1)
      pars_obs_d$X3 <- list(beta = logit(0.6))
    }

    dat.r <- data.frame(U1 = rnorm(n.r,0,1),
                        U2 = rnorm(n.r,0,0.35))
    dat.r <- as.data.table(causl:::rfrugalParam(n = n.r, formulas = forms_exp_d, family = family, pars = pars_exp_d, dat = dat.r))
    dat.r[,  NCO.strong := 1 + X2 + X3 + U1 + rnorm(n.r,0,2)]
    dat.r[,  NCO.weak := 1 + X2 + X3 + U2 + rnorm(n.r,0,2)]

    if (U.diff) {
      dat.o <- data.frame(U1 = rnorm(n.o,U1.mean.shift,1),
                          U2 = rnorm(n.o,U2.mean.shift,0.35))
    } else{
      dat.o <- data.frame(U1 = rnorm(n.o,0,1),
                          U2 = rnorm(n.o,0,0.35))
    }
    dat.o <- as.data.table(causl:::rfrugalParam(n = n.o, formulas = forms_obs_d, family = family, pars = pars_obs_d, dat = dat.o))
    dat.o[,  NCO.strong := 1 + X2 + X3 + U1 + rnorm(n.o,0,2)]
    dat.o[,  NCO.weak := 1 + X2 + X3 + U2 + rnorm(n.o,0,2)]

    dat.stacked <- rbind(dat.r[,S := 1],dat.o[,S := 0])
    dat.stacked[,coeff.U := coeff.U]
    dat.stacked[,seed := seed]
    data.list[[i]] <- as.data.table(dat.stacked[,coeff.U := coeff.U])
  }
  return(data.list)
}

# SIMULATIONS TO MIMIC XI'S
set.seed(124)
B <- 10
n_rct <- 2000
n_rwd <- 2000
est_atmle_glmnet <- matrix(nrow = B, ncol = length(coeff.U.list))
est_atmle_sl <- matrix(nrow = B, ncol = length(coeff.U.list))
est_escvtmle <- matrix(nrow = B, ncol = length(coeff.U.list))
est_rct_only <- matrix(nrow = B, ncol = length(coeff.U.list))

for (b in 1:B) {
  print("run: " %+% b)
  seed.num <- .Random.seed[b]

  # simulate data
  dat.list <- generate.data.list.HTE(n.r = n_rct, n.o = n_rwd,
                                     HTE = TRUE,
                                     seed = seed.num,
                                     coeff.U.list = coeff.U.list,
                                     X.diff = F, U.diff = F)

  for (i in 1:length(coeff.U.list)) {
    data <- as.data.frame(dat.list[[i]])

    # atmle with glmnet
    atmle_glmnet_res <- atmle(data = data,
                              S_node = 10,
                              W_node = c(3, 4, 6),
                              A_node = 5,
                              Y_node = 7,
                              atmle_pooled = TRUE,
                              controls_only = FALSE,
                              theta_method = "glmnet",
                              Pi_method = "glmnet",
                              g_method = "glmnet",
                              theta_tilde_method = "glmnet",
                              Q_method = "glmnet",
                              bias_working_model = "glmnet",
                              pooled_working_model = "glmnet",
                              g_rct = 0.5,
                              family = "gaussian",
                              verbose = FALSE)

    # atmle with super learner
    atmle_sl_res <- atmle(data = data,
                          S_node = 10,
                          W_node = c(3, 4, 6),
                          A_node = 5,
                          Y_node =  7,
                          atmle_pooled = TRUE,
                          controls_only = FALSE,
                          theta_method = "sl3",
                          Pi_method = "sl3",
                          g_method = "sl3",
                          theta_tilde_method = "sl3",
                          Q_method = "glmnet",
                          bias_working_model = "glmnet",
                          pooled_working_model = "glmnet",
                          g_rct = 0.5,
                          family = "gaussian",
                          verbose = FALSE)

    # ES-CVTMLE
    escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
                              data = data,
                              study = "S",
                              covariates = c("X1", "X2", "X3"),
                              treatment_var = "A",
                              treatment = 1,
                              outcome = "Y",
                              pRCT = 0.5,
                              family = "gaussian",
                              Q.SL.library = c("SL.glm"),
                              g.SL.library = c("SL.glm"),
                              Q.discreteSL = TRUE,
                              g.discreteSL = TRUE,
                              V = 10)

    # rct only
    rct_only_res <- rct_only(data = data,
                             S_node = 10,
                             W_node = c(3, 4, 6),
                             A_node = 5,
                             Y_node = 7,
                             g_rct = 0.5,
                             nuisance_method = "sl3",
                             family = "gaussian",
                             verbose = FALSE)

    # collect results
    est_atmle_glmnet[b,i] <- atmle_glmnet_res$est
    est_atmle_sl[b,i] <- atmle_sl_res$est
    est_rct_only[b,i] <- rct_only_res$est
    est_escvtmle[b,i] <- escvtmle_res$ATE$b2v
  }
}

save(list = c("est_atmle_glmnet",
              "est_atmle_sl",
              "est_escvtmle",
              "est_rct_only"), file = "sim_results.RData")
