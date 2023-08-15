library(devtools)

load_all()

`%+%` <- function(a, b) paste0(a, b)

# parameters
bA_truth <- 0.5

fnames <- list.files("out/")

mean(unlist(lapply(fnames, function(fname) {
  obj <- readRDS("out/" %+% fname)
  sapply(obj, function(res) {
    return(ifelse(res$lower <= bA_truth & res$upper >= bA_truth, 1, 0))
  })
})))

