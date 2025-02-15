library(devtools)
library(tmle)
library(torch)
library(hal9001)
library(EScvtmle)
library(furrr)
library(doMC)
library(ggplot2)
load_all()
set.seed(3874) # good
#set.seed(4905090) #bug
plan(multisession, workers = availableCores()-1)
registerDoMC(cores = availableCores()-1)

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
    b <- 0.2+0.1*W1*(1-A)
  } else if (bias == "b") {
    b <- 0.5+3.1*W1*(1-A)+0.8*W3
  }

  # outcome
  Y <- 2.5+0.9*W1+1.1*W2+2.7*W3+ate*A+UY+(1-S)*b

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

# simulate data
data_rct <- sim_data(ate = 1.5,
                     n = 500,
                     rct = TRUE,
                     g_rct = 0.5,
                     bias = "a",
                     controls_only = FALSE)
data_rwd <- sim_data(ate = 1.5,
                     n = 1500,
                     rct = FALSE,
                     g_rct = 0.5,
                     bias = "a",
                     controls_only = FALSE)
data <- rbind(data_rct, data_rwd)
W <- colnames(data)[grep("W", colnames(data))]

# SEQ-A-TMLE
res_seq_atmle <- atmle_torch(data = data,
                             S = "S",
                             W = W,
                             A = "A",
                             Y = "Y",
                             lr = 1e-3,
                             controls_only = FALSE,
                             family = "gaussian",
                             browse = FALSE,
                             parallel = TRUE)
res_seq_atmle_df <- map_dfr(res_seq_atmle, function(.x) {
  data.frame(psi = .x$psi,
             lower = .x$lower,
             upper = .x$upper,
             PnEIC = .x$PnEIC,
             sn = .x$sn)
})
idx <- which.min(res_seq_atmle_df$upper - res_seq_atmle_df$lower)
res_seq_atmle_df$index <- seq_along(res_seq_atmle_df$psi)

# TMLE
res_tmle <- rct_only(data = data,
                     S = "S",
                     W = W,
                     A = "A",
                     Y = "Y",
                     g_rct = 0.5,
                     family = "gaussian",
                     nuisance_method = "glm",
                     verbose = FALSE)
res_seq_atmle_df <- rbind(res_seq_atmle_df,
                          data.frame(psi = res_tmle$est,
                                     lower = res_tmle$lower,
                                     upper = res_tmle$upper,
                                     PnEIC = NA,
                                     sn = NA,
                                     index = length(res_seq_atmle_df$psi) + 1))

ggplot(res_seq_atmle_df, aes(x = index, y = psi)) +
  geom_point(aes(color = ifelse(index == idx, "highlight",
                                ifelse(index == max(index), "tmle", "normal")))) +
  geom_errorbar(aes(ymin = lower, ymax = upper,
                    color = ifelse(index == idx, "highlight",
                                   ifelse(index == max(index), "tmle", "normal"))),
                width = 0.2) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("highlight" = "red", "tmle" = "blue", "normal" = "black"), guide = "none") +
  labs(x = "idx", y = "Psi Estimate",
       title = "Sequence of Psi Estimates with 95% Confidence Intervals") +
  theme_minimal()

abs(res_seq_atmle_df$PnEIC) <= res_seq_atmle_df$sn

