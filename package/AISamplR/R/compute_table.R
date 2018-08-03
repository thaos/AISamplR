
compute_weight_table <- 
  function(loglik_table, denom_table)
  {
    denom_table[is.infinite(denom_table) & denom_table < 0] <-
      min(denom_table[is.finite(denom_table) & denom_table < 0])
    weight_table <- exp(loglik_table - denom_table)
    weight_table[is.infinite(weight_table)] <-
      max(weight_table[is.finite(weight_table)]) + 1
    weight_table <- weight_table / max(weight_table)
    weight_table <- weight_table / sum(weight_table)
    if(median(weight_table) == 0){
      warning("A majority of weights are equal to zero \n
              => could be preferable to modify the sigma parameter of gen_xs_rcpp")
    }
    return(weight_table)
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
    lp_arr <- structure(lp_2dmat, dim = c(T, M, N))
    return(lp_arr)
  }

compute_loglik_table.externalptr  <- 
  function (logposterior, xs_chain, D, T, N, M)
  {
    compute_loglik_table_rcpp(lp = logposterior,
                              x = xs_chain,
                              D = 2, T = T, N = N, M = M)
  }

