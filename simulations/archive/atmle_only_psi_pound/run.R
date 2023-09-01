library(devtools)

load_all()

`%+%` <- function(a, b) paste0(a, b)

# parameters
bA_truth <- 0.5
n_times <- 10
p_rct <- 0.67

estimates <- vector(mode = "list", length = n_times)

for (i in 1:n_times) {
  data <- generate_data(N=500, p_rct=p_rct, bA=bA_truth)
  res <- atmle_tmle(data,
                    S_node = 1,
                    W_node = c(2, 3, 4, 5),
                    A_node = 6,
                    Y_node = 7,
                    target_Pi = TRUE,
                    g_rct=p_rct)
  estimates[[i]] <- res
}

saveRDS(estimates, "out/atmle_only_psi_pound_" %+% gsub(" ", "_", timestamp(prefix = "", suffix = "", quiet = TRUE)) %+% ".RDS")
