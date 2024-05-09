set.seed(123)
n <- 500
S <- rbinom(n, 1, 0.5)
W <- rnorm(n)
A <- rbinom(n, 1, 0.5)
Y <- rnorm(n)
data <- data.frame(S, W, A, Y)

load_all()
res <- atmle(data = data,
             S_node = 1,
             W_node = 2,
             A_node = 3,
             Y_node = 4,
             controls_only = FALSE,
             family = "gaussian")
