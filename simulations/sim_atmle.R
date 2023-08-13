library(parallel)
library(pbmcapply)

times <- seq(1, 200)

final_res <- pbmclapply(times, function(i) {
  data <- generate_data(N=500, p_rct=0.67, bA=0.7)
  res <- atmle(data, g_rct=0.67, max_iter=1, eps=1e-5)

  return(res)
}, mc.cores = 5)

saveRDS(final_res, "atmle_res.RDS")

coverage <- sapply(final_res, function(res) {
  return(ifelse(res$lower <= 0.7 & res$upper >= 0.7, 1, 0))
})
mean(coverage)
