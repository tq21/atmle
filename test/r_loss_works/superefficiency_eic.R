library(ggpubr)
load_all()
source("utils.R")

# parameters -------------------------------------------------------------------
B <- 50
n <- 2000
truth <- 1.5
nuisance_method <- "glm"
working_model <- "lasso"
var_eic_nonparam <- rep(NA, B)
var_eic_atmle_no_bias <- rep(NA, B)
var_eic_atmle_simple_bias <- rep(NA, B)
var_eic_atmle_complex_bias <- rep(NA, B)
var_eic_atmle_misspecify_bias <- rep(NA, B)
var_eic_escvtmle_no_bias <- rep(NA, B)
var_eic_escvtmle_simple_bias <- rep(NA, B)
var_eic_escvtmle_complex_bias <- rep(NA, B)
var_eic_escvtmle_misspecify_bias <- rep(NA, B)

for (i in 1:B) {
  print(i)
  data <- generate_data(n, truth, 0)
  S_node = 1
  W_node = c(2, 3, 4, 5)
  A_node = 6
  Y_node = 7

  nonparam_res <- rct_only(data,
                           S_node,
                           W_node,
                           A_node,
                           Y_node,
                           nuisance_method=nuisance_method,
                           working_model=working_model,
                           p_rct=0.5,
                           verbose=FALSE)
  var_eic_nonparam[i] <- nonparam_res$var_eic

  # no bias
  data <- generate_data(n, truth, 0)
  S_node = 1
  W_node = c(2, 3, 4, 5)
  A_node = 6
  Y_node = 7

  atmle_no_bias <- atmle(data,
                         S_node,
                         W_node,
                         A_node,
                         Y_node,
                         nuisance_method=nuisance_method,
                         working_model=working_model,
                         p_rct=0.5,
                         verbose=FALSE)
  var_eic_atmle_no_bias[i] <- atmle_no_bias$var_eic
  data$S <- 1 - data$S
  escvtmle_no_bias <- ES.cvtmle(txinrwd = TRUE,
                                data = data,
                                study = "S",
                                covariates = c("W1", "W2", "W3", "W4"),
                                treatment_var = "A",
                                treatment = 1,
                                outcome = "Y",
                                pRCT = 0.5,
                                family = "gaussian",
                                Q.SL.library = c("SL.glm"),
                                g.SL.library = c("SL.glm"),
                                Q.discreteSL = TRUE,
                                g.discreteSL = TRUE,
                                V = 5)
  var_eic_escvtmle_no_bias[i] <- as.numeric(((escvtmle_no_bias$CI$b2v[1]-escvtmle_no_bias$ATE$b2v)/1.96)^2*n)

  # simple bias
  data <- generate_data(n, truth, 0.8)
  S_node = 1
  W_node = c(2, 3, 4, 5)
  A_node = 6
  Y_node = 7

  atmle_simple_bias <- atmle(data,
                             S_node,
                             W_node,
                             A_node,
                             Y_node,
                             nuisance_method=nuisance_method,
                             working_model=working_model,
                             p_rct=0.5,
                             verbose=FALSE)
  var_eic_atmle_simple_bias[i] <- atmle_simple_bias$var_eic
  data$S <- 1 - data$S
  escvtmle_simple_bias <- ES.cvtmle(txinrwd = TRUE,
                                    data = data,
                                    study = "S",
                                    covariates = c("W1", "W2", "W3", "W4"),
                                    treatment_var = "A",
                                    treatment = 1,
                                    outcome = "Y",
                                    pRCT = 0.5,
                                    family = "gaussian",
                                    Q.SL.library = c("SL.glm"),
                                    g.SL.library = c("SL.glm"),
                                    Q.discreteSL = TRUE,
                                    g.discreteSL = TRUE,
                                    V = 5)
  var_eic_escvtmle_simple_bias[i] <- as.numeric(((escvtmle_simple_bias$CI$b2v[1]-escvtmle_simple_bias$ATE$b2v)/1.96)^2*n)

  # complex bias
  data <- generate_data(n, truth, "complex")
  S_node = 1
  W_node = c(2, 3, 4, 5)
  A_node = 6
  Y_node = 7

  atmle_complex_bias <- atmle(data,
                              S_node,
                              W_node,
                              A_node,
                              Y_node,
                              nuisance_method=nuisance_method,
                              working_model=working_model,
                              p_rct=0.5,
                              verbose=FALSE)
  var_eic_atmle_complex_bias[i] <- atmle_complex_bias$var_eic
  data$S <- 1 - data$S
  escvtmle_complex_bias <- ES.cvtmle(txinrwd = TRUE,
                                     data = data,
                                     study = "S",
                                     covariates = c("W1", "W2", "W3", "W4"),
                                     treatment_var = "A",
                                     treatment = 1,
                                     outcome = "Y",
                                     pRCT = 0.5,
                                     family = "gaussian",
                                     Q.SL.library = c("SL.glm"),
                                     g.SL.library = c("SL.glm"),
                                     Q.discreteSL = TRUE,
                                     g.discreteSL = TRUE,
                                     V = 5)
  var_eic_escvtmle_complex_bias[i] <- as.numeric(((escvtmle_complex_bias$CI$b2v[1]-escvtmle_complex_bias$ATE$b2v)/1.96)^2*n)

  # misspecified bias
  data <- generate_data(n, truth, "misspecify")
  S_node = 1
  W_node = c(2, 3, 4, 5)
  A_node = 6
  Y_node = 7

  atmle_misspecify_bias <- atmle(data,
                                 S_node,
                                 W_node,
                                 A_node,
                                 Y_node,
                                 nuisance_method=nuisance_method,
                                 working_model=working_model,
                                 p_rct=0.5,
                                 verbose=FALSE)
  var_eic_atmle_misspecify_bias[i] <- atmle_misspecify_bias$var_eic
  data$S <- 1 - data$S
  escvtmle_misspecify_bias <- ES.cvtmle(txinrwd = TRUE,
                                        data = data,
                                        study = "S",
                                        covariates = c("W1", "W2", "W3", "W4"),
                                        treatment_var = "A",
                                        treatment = 1,
                                        outcome = "Y",
                                        pRCT = 0.5,
                                        family = "gaussian",
                                        Q.SL.library = c("SL.glm"),
                                        g.SL.library = c("SL.glm"),
                                        Q.discreteSL = TRUE,
                                        g.discreteSL = TRUE,
                                        V = 5)
  var_eic_escvtmle_misspecify_bias[i] <- as.numeric(((escvtmle_misspecify_bias$CI$b2v[1]-escvtmle_misspecify_bias$ATE$b2v)/1.96)^2*n)
}

data_n <- data.frame(run = 1:B,
                     var_eic_nonparam = var_eic_nonparam,
                     var_eic_atmle_no_bias = var_eic_atmle_no_bias,
                     var_eic_atmle_simple_bias = var_eic_atmle_simple_bias,
                     var_eic_atmle_complex_bias = var_eic_atmle_complex_bias,
                     var_eic_atmle_misspecify_bias = var_eic_atmle_misspecify_bias,
                     var_eic_escvtmle_no_bias = var_eic_escvtmle_no_bias,
                     var_eic_escvtmle_simple_bias = var_eic_escvtmle_simple_bias,
                     var_eic_escvtmle_complex_bias = var_eic_escvtmle_complex_bias,
                     var_eic_escvtmle_misspecify_bias = var_eic_escvtmle_misspecify_bias)

p_atmle <- ggplot(data_n) +
  geom_point(aes(x = run, y = var_eic_nonparam, color = "Nonparametric")) +
  geom_line(aes(x = run, y = var_eic_nonparam, color = "Nonparametric")) +
  geom_point(aes(x = run, y = var_eic_atmle_no_bias, color = "no bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_no_bias, color = "no bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_simple_bias, color = "simple bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_simple_bias, color = "simple bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_complex_bias, color = "complex bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_complex_bias, color = "complex bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_misspecify_bias, color = "misspecified bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_misspecify_bias, color = "misspecified bias")) +
  scale_y_continuous(breaks = seq(0, 30, 10), limits = c(0, 30)) +
  labs(title = "",
       x = "ATMLE",
       y = "Variance of EIC",
       color = "Scenario") +
  scale_color_manual(values = c("Nonparametric" = "black",
                                "no bias" = "blue",
                                "simple bias" = "orange",
                                "complex bias" = "red",
                                "misspecified bias" = "pink")) +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_atmle_focus <- ggplot(data_n) +
  geom_point(aes(x = run, y = var_eic_atmle_no_bias, color = "no bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_no_bias, color = "no bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_simple_bias, color = "simple bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_simple_bias, color = "simple bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_complex_bias, color = "complex bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_complex_bias, color = "complex bias")) +
  scale_y_continuous(breaks = seq(3, 7, 1), limits = c(3, 7)) +
  labs(title = "",
       x = "ATMLE",
       y = "Variance of EIC",
       color = "Scenario") +
  scale_color_manual(values = c("no bias" = "blue",
                                "simple bias" = "orange",
                                "complex bias" = "red")) +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_escvtmle <- ggplot(data_n) +
  geom_point(aes(x = run, y = var_eic_nonparam, color = "Nonparametric")) +
  geom_line(aes(x = run, y = var_eic_nonparam, color = "Nonparametric")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_no_bias, color = "no bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_no_bias, color = "no bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_simple_bias, color = "simple bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_simple_bias, color = "simple bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_complex_bias, color = "complex bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_complex_bias, color = "complex bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_misspecify_bias, color = "misspecified bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_misspecify_bias, color = "misspecified bias")) +
  scale_y_continuous(breaks = seq(0, 30, 10), limits = c(0, 30)) +
  labs(title = "",
       x = "ESCVTMLE",
       y = "Variance of EIC",
       color = "Scenario") +
  scale_color_manual(values = c("Nonparametric" = "black",
                                "no bias" = "blue",
                                "simple bias" = "orange",
                                "complex bias" = "red",
                                "misspecified bias" = "pink")) +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_escvtmle_focus <- ggplot(data_n) +
  geom_point(aes(x = run, y = var_eic_escvtmle_no_bias, color = "no bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_no_bias, color = "no bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_simple_bias, color = "simple bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_simple_bias, color = "simple bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_complex_bias, color = "complex bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_complex_bias, color = "complex bias")) +
  scale_y_continuous(breaks = seq(3, 7, 1), limits = c(3, 7)) +
  labs(title = "",
       x = "ESCVTMLE",
       y = "Variance of EIC",
       color = "Scenario") +
  scale_color_manual(values = c("no bias" = "blue",
                                "simple bias" = "orange",
                                "complex bias" = "red")) +
  theme_minimal() +
  theme(text = element_text(size = 16))

plt_1 <- ggarrange(p_atmle, p_atmle_focus,
                 nrow = 1, ncol = 2, common.legend = TRUE)
plt_2 <- ggarrange(p_atmle_focus, p_escvtmle_focus,
                   nrow = 1, ncol = 2, common.legend = TRUE)
ggsave(filename = "superefficiency_eic_atmle.pdf", plot = plt_1, device = "pdf",
       path = "plot", width = 12, height = 6, dpi = 300)
ggsave(filename = "superefficiency_eic_atmle_escvtmle.pdf", plot = plt_2, device = "pdf",
       path = "plot", width = 12, height = 6, dpi = 300)
