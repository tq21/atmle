#' @title Make repeated measures data for pooled logistic regression
#' for estimating conditional fT_tilde_ilure-time hazard function.
#'
#' @param W A matrix of covariates.
#' @param T_tilde A vector of time-points.
#'
#' @export
make_rep_data <- function(W, T_tilde) {

  # make id
  id <- seq(length(T_tilde))

  # make repeated data for pooled logistic regression
  rep_data <- map_dfr(id, function(i) {
    Wi <- W[i, , drop = FALSE]
    T_tilde_i <- T_tilde[i]
    id <- id[i]
    idi_exp <- rep(id, T_tilde_i)
    T_tilde_i_exp <- seq(T_tilde_i)
    wi_exp <- do.call(rbind, replicate(T_tilde_i, Wi, simplify = FALSE))
    ind <- c(rep(0, T_tilde_i - 1), 1)

    return(data.frame(id = idi_exp, T_tilde_i = T_tilde_i_exp, wi_exp, ind = ind))
  })

  rownames(rep_data) <- NULL

  return(list(id = rep_data$id,
              ind = rep_data$ind,
              data = rep_data[, !names(rep_data) %in% c("id", "ind"), drop = FALSE]))
}
