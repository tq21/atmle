library(ggpubr)
load_all()
source("utils.R")

# parameters -------------------------------------------------------------------
B <- 50
n <- 500
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
  data <- generate_data(500, 1.5, 0)
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
  data <- generate_data(500, 1.5, 0)
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
  data <- generate_data(500, 1.5, 0.8)
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
  data <- generate_data(500, 1.5, "complex")
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
  data <- generate_data(500, 1.5, "misspecify")
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
  geom_point(aes(x = run, y = var_eic_atmle_no_bias, color = "ATMLE no bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_no_bias, color = "ATMLE no bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_simple_bias, color = "ATMLE simple bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_simple_bias, color = "ATMLE simple bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_complex_bias, color = "ATMLE complex bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_complex_bias, color = "ATMLE complex bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_misspecify_bias, color = "ATMLE misspecified bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_misspecify_bias, color = "ATMLE misspecified bias")) +
  scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 40)) +
  labs(title = "",
       x = "ATMLE",
       y = "Variance of EIC",
       color = "Method") +
  scale_color_manual(values = c("Nonparametric" = "black",
                                "ATMLE no bias" = "blue",
                                "ATMLE simple bias" = "orange",
                                "ATMLE complex bias" = "red",
                                "ATMLE misspecified bias" = "pink")) +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_atmle_focus <- ggplot(data_n) +
  geom_point(aes(x = run, y = var_eic_atmle_no_bias, color = "ATMLE no bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_no_bias, color = "ATMLE no bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_simple_bias, color = "ATMLE simple bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_simple_bias, color = "ATMLE simple bias")) +
  geom_point(aes(x = run, y = var_eic_atmle_complex_bias, color = "ATMLE complex bias")) +
  geom_line(aes(x = run, y = var_eic_atmle_complex_bias, color = "ATMLE complex bias")) +
  scale_y_continuous(breaks = seq(3, 8, 1), limits = c(3, 8)) +
  labs(title = "",
       x = "ATMLE",
       y = "Variance of EIC",
       color = "Method") +
  scale_color_manual(values = c("Nonparametric" = "black",
                                "ATMLE no bias" = "blue",
                                "ATMLE simple bias" = "orange",
                                "ATMLE complex bias" = "red")) +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_escvtmle <- ggplot(data_n) +
  geom_point(aes(x = run, y = var_eic_nonparam, color = "Nonparametric")) +
  geom_line(aes(x = run, y = var_eic_nonparam, color = "Nonparametric")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_no_bias, color = "ESCVTMLE no bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_no_bias, color = "ESCVTMLE no bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_simple_bias, color = "ESCVTMLE simple bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_simple_bias, color = "ESCVTMLE simple bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_complex_bias, color = "ESCVTMLE complex bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_complex_bias, color = "ESCVTMLE complex bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_misspecify_bias, color = "ESCVTMLE misspecified bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_misspecify_bias, color = "ESCVTMLE misspecified bias")) +
  scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 40)) +
  labs(title = "",
       x = "ESCVTMLE",
       y = "Variance of EIC",
       color = "Method") +
  scale_color_manual(values = c("Nonparametric" = "black",
                                "ESCVTMLE no bias" = "blue",
                                "ESCVTMLE simple bias" = "orange",
                                "ESCVTMLE complex bias" = "red",
                                "ESCVTMLE misspecified bias" = "pink")) +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_escvtmle_focus <- ggplot(data_n) +
  geom_point(aes(x = run, y = var_eic_escvtmle_no_bias, color = "ESCVTMLE no bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_no_bias, color = "ESCVTMLE no bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_simple_bias, color = "ESCVTMLE simple bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_simple_bias, color = "ESCVTMLE simple bias")) +
  geom_point(aes(x = run, y = var_eic_escvtmle_complex_bias, color = "ESCVTMLE complex bias")) +
  geom_line(aes(x = run, y = var_eic_escvtmle_complex_bias, color = "ESCVTMLE complex bias")) +
  scale_y_continuous(breaks = seq(3, 8, 1), limits = c(3, 8)) +
  labs(title = "",
       x = "ESCVTMLE",
       y = "Variance of EIC",
       color = "Method") +
  scale_color_manual(values = c("ESCVTMLE no bias" = "blue",
                                "ESCVTMLE simple bias" = "orange",
                                "ESCVTMLE complex bias" = "red")) +
  theme_minimal() +
  theme(text = element_text(size = 16))

plt <- ggarrange(p_atmle_focus, p_escvtmle_focus,
                 nrow = 2, ncol = 2, common.legend = TRUE)
# ggsave(filename = "superefficiency_eic.pdf", plot = plt, device = "pdf",
#        path = "plot", width = 12, height = 6, dpi = 300)
