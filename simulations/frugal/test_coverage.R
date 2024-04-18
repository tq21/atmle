# check coverage
B <- 200
ate <- 0.2
atmle_cover <- numeric(B)
for (i in 1:B) {
  data <- rbind(sim_data(2000, 0.5, TRUE), sim_data(1000, 0.5, FALSE))
  S_node <- 1
  W_node <- 2:4
  A_node <- 5
  Y_node <- 6

  atmle_res <- atmle(data = data,
                     S_node = S_node,
                     W_node = W_node,
                     A_node = A_node,
                     Y_node = Y_node,
                     atmle_pooled = TRUE,
                     controls_only = FALSE,
                     theta_method = "glm",
                     Pi_method = "glm",
                     g_method = "glm",
                     theta_tilde_method = "glm",
                     Q_method = "glm",
                     bias_working_model = "glmnet",
                     pooled_working_model = "HAL",
                     g_rct = 0.5,
                     family = "gaussian",
                     verbose = FALSE)

  atmle_cover[i] <- (atmle_res$lower <= ate) & (atmle_res$upper >= ate)

  if (atmle_cover[i] == 1) {
    print("covered")
  } else {
    print("not covered")
  }
}
mean(atmle_cover)


print("atmle coverage: " %+% mean(atmle_cover))
print("escvtmle coverage: " %+% mean(escvtmle_cover))
print("atmle mse: " %+% (round((mean(atmle_est)-ate)^2+var(atmle_est),5)))
print("escvtmle mse: " %+% (round((mean(escvtmle_est)-ate)^2+var(escvtmle_est),5)))
