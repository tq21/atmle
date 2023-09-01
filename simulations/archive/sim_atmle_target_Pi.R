library(parallel)
library(pbmcapply)
library(devtools)

load_all()

`%+%` <- function(a, b) paste0(a, b)

# parameters
bA_truth <- 0.8
n_times <- 500
p_rct <- 0.67
mc.cores <- 20

estimates <- pbmclapply(seq(1, n_times), function(i) {
  data <- generate_data(N=500, p_rct=p_rct, bA=bA_truth)
  res <- atmle(data,
               S_node = 1,
               W_node = c(2, 3, 4, 5),
               A_node = 6,
               Y_node = 7,
               target_Pi = TRUE,
               g_rct=p_rct)
  return(res)
}, mc.cores = mc.cores)

# compute coverage of 95% CIs --------------------------------------------------
coverage <- sapply(estimates, function(res) {
  if (length(res) == 3) {
    return(ifelse(res$lower <= bA_truth & res$upper >= bA_truth, 1, 0))
  } else {
    return(NA)
  }
})
coverage <- mean(coverage, na.rm = TRUE)

# compute bias -----------------------------------------------------------------
bias <- mean(sapply(estimates, function(res) {
  if (length(res) == 3) {
    return(res$est - bA_truth)
  } else {
    return(NA)
  }
}), na.rm = TRUE)

# compute variance -------------------------------------------------------------
vars <- var(sapply(estimates, function(res) {
  if (length(res) == 3) {
    return(res$est)
  } else {
    return(NA)
  }
}), na.rm = TRUE)

# compute MSE
mse <- bias^2 + vars

# save output
out_obj <- list(estimates = estimates,
                coverage = coverage,
                bias = bias,
                vars = vars,
                mse = mse)
saveRDS(out_obj, "out/atmle_target_Pi_" %+% gsub(" ", "_", timestamp(prefix = "", suffix = "", quiet = TRUE)) %+% ".RDS")
