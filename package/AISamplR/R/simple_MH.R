simple_MH <- function(logposterior, d, T = 100, mu = NULL, sigma){
  if(d != length(sigma)) stop("d != length(sigma)")
  x <- matrix(0, nrow = T, ncol=d)
  if(!is.null(mu)) x[1, ] <- mu
  proposition  <- rnorm(d, mean = x[1, ], sd = sigma)
  if(d != length(proposition)) stop("d != length(proposition)")
  prev_lposterior <- logposterior(x[1, ])
  for (i in 2:T){
    proposition  <- rnorm(d, mean = x[i-1, ], sd = sigma)
    #acceptance prob
    cur_lposterior <- logposterior(proposition)
    rho <- exp(cur_lposterior - prev_lposterior)
    alpha <- min(c(1,rho))
    # MH TEST
    u <- runif(1)
    if (u <= alpha){
      x[i, ] <- proposition
      prev_lposterior <- cur_lposterior
    } else{
      x[i, ] <- x[i-1, ]
    }
  }
  return(x)
}
