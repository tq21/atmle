library(parallel)
library(pbmcapply)

times <- seq(1, 200)

final_res <- pbmclapply(times, function(i) {
  data <- generate_data(N=500, p_rct=0.67, bA=0.7)
  res <- parametric(data, g_rct = 0.67)

  return(res)
}, mc.cores = 5)



saveRDS(final_res, "result.RDS")

coverage <- sapply(final_res, function(res) {
  return(ifelse(res$lower <= 0.7 & res$upper >= 0.7, 1, 0))
})
mean(coverage)

res_vec <- vector(length = 200)

for (i in 1:200) {
  res <- final_res[[i]]
  if (!is.null(names(res))) {
    res_vec[i] <- ifelse(res$lower <= 0.7 & res$upper >= 0.7, 1, 0)
  } else {
    res_vec[i] <- 1
  }
}
