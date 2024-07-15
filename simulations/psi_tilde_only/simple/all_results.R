source("utils.R")
source("utils_plot.R")

set.seed(123)

n_seq <- seq(500, 2000, 500)
t0 <- 3

# 001
source("sim_data_001.R")
res_list <- readRDS("out/res_list_001.RDS")
plots_001 <- plt_fun()

# 002
source("sim_data_002.R")
res_list <- readRDS("out/res_list_002.RDS")
plots_002 <- plt_fun()

# 003
source("sim_data_003.R")
res_list <- readRDS("out/res_list_003.RDS")
plots_003 <- plt_fun()

# 004
source("sim_data_004.R")
res_list <- readRDS("out/res_list_004.RDS")
plots_004 <- plt_fun()

plt_list <- list(plots_001, plots_002, plots_003, plots_004)
plt_list <- unlist(plt_list, recursive = FALSE)

plt <- ggarrange(plotlist = plt_list, nrow = 4, ncol = 4, common.legend = TRUE)
ggsave(filename = "all.pdf", plot = plt, device = "pdf",
       path = "figs", width = 16, height = 12, dpi = 300)
