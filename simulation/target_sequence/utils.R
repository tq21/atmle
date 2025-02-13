.libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.2",
            .libPaths()))
library(devtools)
library(tmle)
library(torch)
library(hal9001)
library(EScvtmle)
load_all()
source("sim_data.R")
set.seed(123)
run <- function(B,
                bias,
                g_rct,
                controls_only,
                ate,
                n_rct_seq,
                n_rwd_seq) {

  res_df <- map_dfr(seq(length(n_rct_seq)), function(j) {
    map_dfr(seq(B), function(b) {
      cat("n =", n_rct_seq[j] + n_rwd_seq[j], "b =", b, "\n")

      # simulate data
      data_rct <- sim_data(ate = ate,
                           n = n_rct_seq[j],
                           rct = TRUE,
                           g_rct = g_rct,
                           bias = bias,
                           controls_only = controls_only)
      data_rwd <- sim_data(ate = ate,
                           n = n_rwd_seq[j],
                           rct = FALSE,
                           g_rct = g_rct,
                           bias = bias,
                           controls_only = controls_only)
      data <- rbind(data_rct, data_rwd)
      W <- colnames(data)[grep("W", colnames(data))]
      res_df_tmp <- data.frame(bias = bias,
                               n = n_rct_seq[j] + n_rwd_seq[j],
                               b = b,
                               est_name = numeric(3),
                               psi = numeric(3),
                               lower = numeric(3),
                               upper = numeric(3))

      # SEQ-A-TMLE
      res_seq_atmle <- atmle_torch(data = data,
                                   S = "S",
                                   W = W,
                                   A = "A",
                                   Y = "Y",
                                   controls_only = controls_only,
                                   family = "gaussian",
                                   browse = FALSE)
      res_seq_atmle_df <- map_dfr(res_seq_atmle, function(.x) {
        data.frame(psi = .x$psi,
                   lower = .x$lower,
                   upper = .x$upper)
      })
      idx <- which.min(res_seq_atmle_df$upper - res_seq_atmle_df$lower)
      res_df_tmp[1, 4:7] <- c("A-TMLE", res_seq_atmle_df[1, ])
      res_df_tmp[2, 4:7] <- c("SEQ-A-TMLE", res_seq_atmle_df[idx, ])

      # ES-CVTMLE
      res_escvtmle <- ES.cvtmle(txinrwd = !controls_only,
                                data = data,
                                study = "S",
                                covariates = W,
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
      res_df_tmp[3, 4:7] <- c("ES-CVTMLE",
                              res_escvtmle$ATE$b2v,
                              as.numeric(res_escvtmle$CI$b2v[1]),
                              as.numeric(res_escvtmle$CI$b2v[2]))

      return(res_df_tmp)
    })
  })

  return(res_df)
}
