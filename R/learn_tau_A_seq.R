learn_tau_A_seq <- function(S,
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
                            parallel) {

  # use HAL to learn tau_A
  pseudo_outcome <- (Y-theta)/(A-g1W)
  pseudo_weights <- (A-g1W)^2*weights
  basis_list <- enumerate_basis(
    x = as.matrix(W),
    max_degree = enumerate_basis_args$max_degree,
    smoothness_orders = enumerate_basis_args$smoothness_orders
  )
  hal_design <- make_design_matrix(X = as.matrix(W), blist = basis_list)

  # run optimization routine to target betas in each working model
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
    loss_fn = r_loss_cate,
    lr = lr,
    max_iter = max_iter,
    verbose = verbose,
    device = device,
    tolerance = tolerance,
    patience = patience,
    parallel = parallel
  )

  return(list(beta_list = beta_list,
              hal_design = hal_design))
}
