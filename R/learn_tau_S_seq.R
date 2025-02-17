learn_tau_S_seq <- function(S,
                            W,
                            A,
                            Y,
                            theta,
                            Q_bar,
                            g1W,
                            Pi1WA,
                            weights,
                            enumerate_basis_args,
                            lr,
                            max_iter,
                            verbose,
                            device,
                            tolerance,
                            patience,
                            parallel,
                            Pi_iter_target) {

  # use HAL to learn tau_S
  pseudo_outcome <- (Y-Q_bar)/(A-Pi1WA)
  pseudo_weights <- (S-Pi1WA)^2*weights
  basis_list <- enumerate_basis(
    x = as.matrix(cbind(W, A = A)),
    max_degree = enumerate_basis_args$max_degree,
    smoothness_orders = enumerate_basis_args$smoothness_orders
  )
  hal_design <- make_design_matrix(X = as.matrix(cbind(W, A = A)), blist = basis_list)
  hal_design_A1 <- make_design_matrix(X = as.matrix(cbind(W, A = 1)), blist = basis_list)
  hal_design_A0 <- make_design_matrix(X = as.matrix(cbind(W, A = 0)), blist = basis_list)

  # run optimization routine to target betas in each working model
  if (Pi_iter_target) {
    beta_list <- torch_optim_routine_with_Pi(
      S = S,
      A = A,
      Y = Y,
      theta = theta,
      Q_bar = Q_bar,
      g1W = g1W,
      Pi1WA = Pi1WA,
      hal_design = hal_design,
      pseudo_outcome = pseudo_outcome,
      pseudo_weights = pseudo_weights,
      loss_fn = r_loss_care_with_Pi,
      lr = lr,
      max_iter = max_iter,
      verbose = verbose,
      device = device,
      tolerance = tolerance,
      patience = patience,
      parallel = parallel
    )
  } else {
    beta_list <- torch_optim_routine(
      S = S,
      A = A,
      Y = Y,
      theta = theta,
      Q_bar = Q_bar,
      g1W = g1W,
      Pi1WA = Pi1WA,
      hal_design = hal_design,
      pseudo_outcome = pseudo_outcome,
      pseudo_weights = pseudo_weights,
      loss_fn = r_loss_care,
      lr = lr,
      max_iter = max_iter,
      verbose = verbose,
      device = device,
      tolerance = tolerance,
      patience = patience,
      parallel = parallel
    )
  }

  return(list(beta_list = beta_list,
              hal_design = hal_design,
              hal_design_A1 = hal_design_A1,
              hal_design_A0 = hal_design_A0))
}
