lambda_tmle <- function(A,
                        t0) {

  unique_T_tilde <- sort(unique(T_tilde))
  t0_idx <- which(t0 == unique_T_tilde)
  t_idx_grid <- seq(t0_idx)

  # clever covariate
  get_H <- function(a) {

    if (a == 1) {
      lambda_pred <- lambda$A1
      G_bar_pred <- G_bar$A1
    } else if (a == 0) {
      lambda_pred <- lambda$A0
      G_bar_pred <- G_bar$A0
    }

    surv_t0 <- map_dbl(seq(length(A)), function(.id) {
      hazards <- lambda_pred[.id == lambda$id][seq(t0_idx)]
      if (length(hazards) == 1) {
        return(hazards)
      } else {
        return(prod(1 - hazards[1:(length(hazards) - 1)])*hazards[length(hazards)])
      }
    })

    surv_t <- map(t_idx_grid, function(.t_idx) {
      map_dbl(seq(length(A)), function(.id) {
        hazards <- lambda_pred[.id == lambda$id][seq(.t_idx)]
        if (length(hazards) == 1) {
          return(hazards)
        } else {
          return(prod(1 - hazards[1:(length(hazards) - 1)])*hazards[length(hazards)])
        }
      })
    })

    H <- unlist(map(t_idx_grid, function(.t_idx) {
      surv_t_cur <- surv_t[[.t_idx]]
      data.frame(H = -as.numeric(A == a)/(g*G_bar_pred)*(surv_t0/surv_t_cur))
    }))

    return(H)
  }

  # clever covariates
  H_A1 <- get_H(a = 1)
  H_A0 <- get_H(a = 0)
  H_all <- data.frame(H_A1 = H_A1,
                      H_A0 = H_A0,
                      t = rep(unique_T_tilde[t_idx_grid], each = length(A)))

  H_tilde <- g*(1-g)*G_bar$integrate_A*(H_A1-H_A0) # TODO: multiply by phi(W)

  # TODO: TMLE: local least favorable submodel update
}
