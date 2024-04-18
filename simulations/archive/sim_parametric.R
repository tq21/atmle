library(parallel)
library(pbmcapply)
library(devtools)

load_all()

`%+%` <- function(a, b) paste0(a, b)

# parameters
bA_truth <- 0.7
n_times <- 200
p_rct <- 0.67
mc.cores <- 8

estimates <- pbmclapply(seq(1, n_times), function(i) {
  data <- generate_data(N=500, p_rct=p_rct, bA=bA_truth)
  res <- parametric(data,
               S_node = 1,
               W_node = c(2, 3, 4, 5),
               A_node = 6,
               Y_node = 7,
               g_rct=p_rct)
  return(res)
}, mc.cores = mc.cores)

coverage <- sapply(estimates, function(res) {
  return(ifelse(res$lower <= bA_truth & res$upper >= bA_truth, 1, 0))
})
coverage <- mean(coverage, na.rm = TRUE)

# save output
out_obj <- list(estimates = estimates,
                coverage = coverage)
saveRDS(out_obj, "out/atmle_" %+% gsub(" ", "_", timestamp(prefix = "", suffix = "", quiet = TRUE)) %+% ".RDS")
