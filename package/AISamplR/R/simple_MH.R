simple_MH <- function(logposterior, mu, sigma, T = 100){
  UseMethod("simple_MH", logposterior)
}

simple_MH.function <- function(logposterior, mu, sigma, T = 100){
  if(length(mu) != length(sigma)) stop("length(mu) != length(sigma)")
  D <- length(mu)
  x <- matrix(0, nrow = D, ncol = T)
  x[, 1] <- mu
  proposition  <- rnorm(D, mean = x[, 1], sd = sigma)
  if(length(mu) != length(proposition)) stop("length(mu) != length(proposition)")
  prev_lposterior <- logposterior(x[, 1])
  for (i in 2:T){
    proposition  <- rnorm(D, mean = x[, i-1], sd = sigma)
    #acceptance prob
    cur_lposterior <- logposterior(proposition)
    rho <- exp(cur_lposterior - prev_lposterior)
    alpha <- min(c(1,rho))
    # MH TEST
    u <- runif(1)
    if (u <= alpha){
      x[, i] <- proposition
      prev_lposterior <- cur_lposterior
    } else{
      x[, i] <- x[, i-1]
    }
  }
  return(x)
}

simple_MH.externalptr <- function(logposterior, mu, sigma, T = 100){
  simple_MH_rcpp(
    lp = logposterior,
    mu = mu,
    sigma = sigma,
    T = T
  )
}

