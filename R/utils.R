normalize_weights <- function(weights){
  if(all(weights == 0)){
    weights <- rep(1, length(weights))
    warning("all weights are zero, replace by equiprobable weights")
  }
  weights / sum(weights)
}

compute_pitable <- function(x_arr, logposterior){
  apply(x_arr, c(1, 3, 4), logposterior)
}

compute_weights <- function(pi_arr, denom_arr){
  exp(pi_arr - denom_arr)
}

xarr4d_tomatrix <- function(x_arr){
  t(matrix(aperm(x_arr, c(2, 1, 3, 4)), nrow = dim(x_arr)[2]))
}

compute_expectation <- function(x_arr, w_arr, f = identity, ...){
  xmat <- xarr4d_tomatrix(x_arr)
  wvect <- as.vector(w_arr)
  l_fx <- apply(xmat, 1, f, ...)
  l_w <- normalize_weights(wvect)
  print(summary(l_w))
  apply(l_fx, 1, weighted.mean, w = l_w)
}

