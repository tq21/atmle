## code to prepare `simulData` dataset goes here

set.seed(123)

# Log linear model
'outcome_mod_lognorm' <- function(X, A, sigma = 0.25) {
  n <- length(A)
  return( -100 + 27.4 * X[, 3]*A + 13.7 * X[, 4] + 10 * X[, 4] * A + 13.7 * X[, 5] - 10 * X[, 5] * A   + exp(sigma * rnorm(n)))
}

'sample_prop' <- function(X, beta = 2) {
  exp(- beta + 2 * X[, 1] + 0.3 * X[, 2] - 0.4 * X[, 3])
}

'trt_prop' <- function(X) {
  expit(- X[, 1] + 0.4 * X[, 2] - 0.25 * X[, 3] - 0.1 * X[, 4] + 0.1 * X[, 5])
}


'gen_Z' <- function(X) {
  z1 <- exp(X[, 1] / 10)
  z2 <- (X[, 3] + X[, 5] + 20)^2
  z3 <- X[, 2] / (2 + 0.5 * exp(X[, 4]))
  z4 <- (X[, 1] + X[, 4] + 20)^2
  z5 <- 0.5 * X[, 2]*X[, 3] + 0.5 * X[, 5]
  return(cbind(z1, z2, z3, z4, z5))
}

p <- 5
n.p <- 2000
n <- 2e4 ## RCT population size
n2 <-2e5 ## RWE population size

## Generating data
X <- matrix(rnorm(n * p, 1, 1), ncol = p)
X <- scale(X) + 1
Z <- gen_Z(X)
Z <- scale(Z) + 1

eS <- sample_prop(X, beta = 7.7)
eS[eS >= 1] = 0.99
S <- sapply(eS, rbinom, n = 1, size = 1)
S.ind <- which(S == 1)

##--------------------------------------------------------------------
## RCT data
n.trial <- length(S.ind)
X.trial <- X[S.ind, ]
Z.trial <- Z[S.ind, ]
A.trial <- rbinom(n.trial, 1, 0.5) # randomly assign treatment to trial participant
Y.trial <- outcome_mod_lognorm(Z.trial, A.trial) # Incorrect outcome model


##--------------------------------------------------------------------
## Real world data
X <- matrix(rnorm(n2 * p, 1, 0.5), ncol = p)
X <- scale(X) + 1
Z <- gen_Z(X)
Z <- scale(Z) + 1

P.ind <- sample(1:n2, size = n.p) ## RWD id
X.rwe <- X[P.ind, ]
Z.rwe <- Z[P.ind, ]
eA <- trt_prop(X.rwe) # treatment propensity score
A.rwe <- sapply(eA, rbinom, n = 1, size = 1)
Y.rwe <- outcome_mod_lognorm(Z.rwe, A.rwe)


data.trial <- data.frame(Y = Y.trial, A = A.trial, X.trial)
data.rwe <- data.frame(Y = Y.rwe, A = A.rwe, X.rwe)

simulData <- rbind(data.trial, data.rwe)
simulData$delta <- c(rep(1, nrow(X.trial)), rep(0, nrow(X.rwe)))


var_method <- "ic"
theta_method <- "sl3"
Pi_method <- "sl3"
g_method <- "sl3"
theta_tilde_method <- "sl3"
Q_method <- "sl3"
bias_working_model <- "glmnet"
pooled_working_model <- "glmnet"
verbose <- TRUE
family <- "gaussian"
atmle_res <- atmle(data = data,
                   S_node = 4,
                   W_node = c(2, 3),
                   A_node = 1,
                   Y_node = 5,
                   atmle_pooled = TRUE,
                   controls_only = FALSE,
                   theta_method = theta_method,
                   Pi_method = Pi_method,
                   g_method = g_method,
                   theta_tilde_method = theta_tilde_method,
                   Q_method = Q_method,
                   bias_working_model = bias_working_model,
                   pooled_working_model = pooled_working_model,
                   g_rct = 0.5,
                   family = family,
                   verbose = FALSE)

escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
                          data = data,
                          study = "S",
                          covariates = c("X", "V2"),
                          treatment_var = "A",
                          treatment = 1,
                          outcome = "Y",
                          pRCT = 0.5,
                          family = family,
                          Q.SL.library = c("SL.glm"),
                          g.SL.library = c("SL.glm"),
                          Q.discreteSL = TRUE,
                          g.discreteSL = TRUE,
                          V = 5)


IntegrativeHTEcf::IntHTEcf( A, X, X.hte, X.cf, Y, S, nboots=50)
