bias_seq <- seq(0, 2, 0.1)
ate <- 0.2

# get bias as multiples of sd of estimates
(apply(est_rwd, 2, mean)-ate)/apply(est_rwd, 2, sd)

# get mse
rmse_atmle_glmnet <- apply(est_atmle_glmnet, 2, function(x) sqrt(mean((x-ate)^2+var(x))))
rmse_escvtmle <- apply(est_escvtmle, 2, function(x) sqrt(mean((x-ate)^2+var(x))))
rmse_rct_only <- apply(est_rct_only, 2, function(x) sqrt(mean((x-ate)^2+var(x))))

# get bias
bias_atmle_glmnet <- apply(est_atmle_glmnet, 2, mean)-ate
bias_escvtmle <- apply(est_escvtmle, 2, mean)-ate
bias_rct_only <- apply(est_rct_only, 2, mean)-ate

# get variance
var_atmle_glmnet <- apply(est_atmle_glmnet, 2, var)
var_escvtmle <- apply(est_escvtmle, 2, var)
var_rct_only <- apply(est_rct_only, 2, var)

# plot
par(mfrow = c(1, 3))

# 1. relative RMSE
plot(bias_seq, rmse_atmle_glmnet/rmse_rct_only, type = "l", ylab = "relative RMSE", col = "red", ylim = c(0.5, 1.5))
lines(bias_seq, rmse_escvtmle/rmse_rct_only, type = "l", col = "orange")

# 2. bias
plot(bias_seq, bias_atmle_glmnet^2, type = "l", ylab = "bias", col = "red", ylim = c(0, 0.01))
lines(bias_seq, bias_escvtmle^2, type = "l", col = "orange")
lines(bias_seq, bias_rct_only^2, type = "l", col = "black")

# 3. variance
plot(bias_seq, var_atmle_glmnet, type = "l", ylab = "variance", col = "red", ylim = c(0, 0.06))
lines(bias_seq, var_escvtmle, type = "l", col = "orange")
lines(bias_seq, var_rct_only, type = "l", col = "black")
