torch_optim_routine_with_Pi <- function(S,
                                        A,
                                        Y,
                                        theta,
                                        Q_bar,
                                        g1W,
                                        Pi1WA,
                                        hal_design,
                                        pseudo_outcome,
                                        pseudo_weights,
                                        loss_fn,
                                        lr,
                                        max_iter,
                                        verbose,
                                        device,
                                        tolerance,
                                        patience,
                                        parallel) {

  # HAL fit
  fit <- cv.glmnet(x = hal_design,
                   y = pseudo_outcome,
                   weights = pseudo_weights,
                   alpha = 1,
                   nlambda = 50,
                   family = "gaussian",
                   parallel = parallel)
  cv_lambda <- fit$lambda.min
  lambda_seq <- fit$lambda
  lambda_seq <- lambda_seq[lambda_seq <= cv_lambda]
  cv_coefs <- as.numeric(coef(fit, s = "lambda.min"))
  cv_idx <- which(cv_coefs[-1] != 0) # exclude intercept (used for subsetting design matrix)
  cv_coefs <- cv_coefs[cv_coefs != 0] # include intercept (used for initial)
  working_model_seq <- map(lambda_seq, function(lambda) {
    cur_idx <- which(as.numeric(coef(fit, s = lambda))[-1] != 0)
    return(list(lambda = lambda,
                idx = sort(union(cur_idx, cv_idx))))
  })

  beta_list <- future_map(working_model_seq, function(.wm) {
    Pi1WA <- torch_tensor(Pi1WA, device = device)
    intercept <- torch_ones(length(Y), 1, device = device)
    S <- torch_tensor(S, device = device)
    A <- torch_tensor(A, device = device)
    Y <- torch_tensor(Y, device = device)

    # make design matrix for the current working model
    phi <- torch_tensor(as.matrix(hal_design[, .wm$idx, drop = FALSE]),
                        device = device)
    phi <- torch_cat(list(intercept, phi), dim = 2L)
    which_in_cv_wm <- which(.wm$idx %in% cv_idx)

    # beta from CV fit as initial values
    beta_padded <- numeric(length(.wm$idx))
    beta_padded[which_in_cv_wm] <- cv_coefs[-1]
    beta_padded <- c(cv_coefs[1], beta_padded)
    beta <- torch_tensor(beta_padded, requires_grad = TRUE)

    # epsilon to target Pi
    epsilon <- torch_zeros(1, device = device, requires_grad = TRUE)

    # initialize optimizer and LR scheduler
    beta_optim <- optim_adam(list(beta), lr = lr)
    epsilon_optim <- optim_adam(list(epsilon), lr = lr / 10)
    train_losses <- numeric(max_iter)
    no_improve_counter <- 0
    best_loss <- Inf
    beta_scheduler <- lr_reduce_on_plateau(beta_optim, mode = 'min',
                                           patience = patience,
                                           factor = 0.5, verbose = verbose)
    epsilon_scheduler <- lr_reduce_on_plateau(epsilon_optim, mode = 'min',
                                              patience = patience,
                                              factor = 0.5, verbose = verbose)

    ones <- torch_ones(length(S), device = device)
    eps <- torch_tensor(1e-6, device = device)

    theta <- torch_tensor(theta, device = device)
    Q_bar <- torch_tensor(Q_bar, device = device)
    g1W <- torch_tensor(g1W, device = device)

    # optimization
    for (j in seq(max_iter)) {
      beta_optim$zero_grad()
      train_loss <- loss_fn(
        S = S,
        A = A,
        Y = Y,
        theta = theta,
        Q_bar = Q_bar,
        g1W = g1W,
        Pi1WA = Pi1WA,
        phi = phi,
        beta = beta,
        epsilon = epsilon,
        ones = ones,
        eps = eps
      )
      train_loss$backward()
      beta_optim$step()

      epsilon_optim$zero_grad()
      train_loss <- loss_fn(
        S = S,
        A = A,
        Y = Y,
        theta = theta,
        Q_bar = Q_bar,
        g1W = g1W,
        Pi1WA = Pi1WA,
        phi = phi,
        beta = beta,
        epsilon = epsilon,
        ones = ones,
        eps = eps
      )
      train_loss$backward()
      epsilon_optim$step()

      #cat("loss", train_loss$item(), "\n")
      train_losses[j] <- train_loss$item()

      if (train_loss$item() < best_loss - tolerance) {
        best_loss <- train_loss$item()
        no_improve_counter <- 0
      } else {
        no_improve_counter <- no_improve_counter + 1
      }

      if (no_improve_counter >= patience) {
        if (verbose) {
          cat("Early stopping at iteration", j, "with train loss:",
              train_loss$item(), "\n")
        }
        train_losses <- train_losses[1:j]
        break
      }

      if (verbose && j %% 100 == 0) {
        cat("iteration", j, "train loss:", train_loss$item(), "\n")
      }

      beta_scheduler$step(train_loss)
      epsilon_scheduler$step(train_loss)
    }

    return(list(lambda = .wm$lambda,
                idx = .wm$idx,
                beta = as.numeric(beta),
                epsilon = as.numeric(epsilon)))
  }, .progress = TRUE)

  return(beta_list)
}
