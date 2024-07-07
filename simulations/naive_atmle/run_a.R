library(devtools)
library(tmle)
load_all()
source("sim_data.R")
source("naive_atmle.R")
`%+%` <- function(a, b) paste0(a, b)

set.seed(2389745)
B <- 500
n_rct_seq <- c(400, 500, 600, 700, 800)
n_rwd_seq <- 3*n_rct_seq
ate <- 1.5
bias <- "a"
g_rct <- 0.67
all_tmle <- vector(mode = "list", length = length(n_rct_seq))
all_naive_atmle <- vector(mode = "list", length = length(n_rct_seq))
all_atmle <- vector(mode = "list", length = length(n_rct_seq))

for (i in 1:length(n_rct_seq)) {
  n_rct <- n_rct_seq[i]; n_rwd <- n_rwd_seq[i]
  print("RCT n = " %+% n_rct %+% ", RWD n = " %+% n_rwd)
  cur_tmle <- numeric(B)
  cur_naive_atmle <- numeric(B)
  cur_atmle <- numeric(B)

  for (b in 1:B) {
    print(b)
    data_rct <- sim_data(ate = ate,
                         n = n_rct,
                         rct = TRUE,
                         g_rct = g_rct,
                         bias = bias,
                         controls_only = FALSE)
    data_rwd <- sim_data(ate = ate,
                         n = n_rwd,
                         rct = FALSE,
                         g_rct = g_rct,
                         bias = bias,
                         controls_only = FALSE)
    data <- rbind(data_rct, data_rwd)
    res_naive <- naive_atmle(data = data,
                             S = "S",
                             W = c("W1", "W2", "W3"),
                             A = "A",
                             Y = "Y",
                             family = "gaussian",
                             g_method = "glm",
                             theta_tilde_method = "glm")
    res <- atmle(data = data,
                 S = "S",
                 W = c("W1", "W2", "W3"),
                 A = "A",
                 Y = "Y",
                 controls_only = FALSE,
                 family = "gaussian",
                 theta_method = "glm",
                 Pi_method = "glm",
                 g_method = "glm",
                 theta_tilde_method = "glm",
                 bias_working_model = "glmnet",
                 pooled_working_model = "glmnet",
                 verbose = FALSE)

    res_tmle <- nonparametric(data = data,
                              S = "S",
                              W = c("W1", "W2", "W3"),
                              A = "A",
                              Y = "Y",
                              controls_only = FALSE,
                              family = "gaussian",
                              g_method = "glm",
                              Pi_method = "glm",
                              Q_method = "glm",
                              verbose = FALSE)$est
    cur_atmle[b] <- res$est
    cur_naive_atmle[b] <- res_naive
    cur_tmle[b] <- res_tmle
  }

  all_atmle[[i]] <- cur_atmle
  all_naive_atmle[[i]] <- cur_naive_atmle
  all_tmle[[i]] <- cur_tmle
}

save(list = c("all_atmle", "all_naive_atmle", "all_tmle"), file = "out/bias_a.RData")
