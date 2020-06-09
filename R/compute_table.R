
compute_weight_table <- 
  function(loglik_table, logdenom_table)
  {
    if(any(is.infinite(logdenom_table))){
      logdenom_table[is.infinite(logdenom_table) & logdenom_table < 0] <-
        min(logdenom_table[is.finite(logdenom_table) & logdenom_table < 0])
      warning("some negative infinite logdenominators: \n
              replacing by minimmum finite denoninator")
    }
    weight_table <- exp(loglik_table - logdenom_table)
    if(any(is.infinite(weight_table))){
      weight_table[is.infinite(weight_table)] <-
        max(weight_table[is.finite(weight_table)]) + 1
      warning("some infinite weights: \n
              replacing by maximum finite weight + 1")
    }
    return(weight_table)
  }

normalize_weights <- function(weights){
  if(all(weights == 0)){
    weights <- rep(1, length(weights))
    warning("all weights are zero, replaced by equiprobable weights")
  }
  weights / sum(weights)
  weights <- weights / max(weights)
  weights <- weights / sum(weights)
  if(median(weights) == 0){
    warning("A majority of weights are equal to zero \n
            => could be preferable to modify
            the sigma parameter of the proposal")
  }
  return(weights)
}
compute_loglik_table <-
  function(logposterior, xs_chain, D, T, N, M)
  {
    UseMethod("compute_loglik_table", logposterior)
  }

compute_loglik_table.function <- 
  function(logposterior, xs_chain, D, T, N, M)
  {
    xs_2dmat <- matrix(xs_chain, nrow = D) 
    lp_2dmat <- eval_logposterior(xs_2dmat, logposterior) 
    lp_arr <- structure(lp_2dmat, dim = c(T, N, M))
    return(lp_arr)
  }

compute_loglik_table.externalptr  <- 
  function (logposterior, xs_chain, D, T, N, M)
  {
    compute_loglik_table_rcpp(lp = logposterior,
                              x = xs_chain,
                              D = D, T = T, N = N, M = M)
  }

