source("utils.R")
source("utils_plot.R")
source("sim_data_001.R")
res_list <- readRDS("out/res_list_001.RDS")
set.seed(123)

n_seq <- seq(500, 2000, 500)
t0 <- 3
plots <- plt_fun()
plt <- ggarrange(plotlist = plots, nrow = 1, ncol = 6, common.legend = TRUE)
ggsave(filename = "simple_001.pdf", plot = plt, device = "pdf",
       path = "figs", width = 18, height = 4, dpi = 300)
