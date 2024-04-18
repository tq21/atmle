library(parallel)
library(pbmcapply)
library(devtools)

load_all()

`%+%` <- function(a, b) paste0(a, b)

# parameters
bA_truth <- 0.8
n_times <- 30
p_rct <- 0.67


fnames <- list.files("out/")


mean(unlist(lapply(fnames, function(fname) {
  obj <- readRDS("out/" %+% fname)
  sapply(obj$estimates, function(res) {
    return(ifelse(res$lower <= bA_truth & res$upper >= bA_truth, 1, 0))
  })
})))

