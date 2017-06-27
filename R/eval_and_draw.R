eval_proposal <- function(x, mu, sigma){
  ans <- exp(sum(dnorm(x = x, mean = mu, sd = sigma, log = TRUE)))# + .Machine$double.xmin
  #if(any(is.na(ans) | ans == 0)) browser()
  ans
}

eval_proposals <- function(xmat, mu, sigma){
  apply(xmat, 2, eval_proposal, mu = mu, sigma = sigma)
}

eval_logposterior <- function(xmat, logposterior){
  apply(xmat, 2, logposterior)
}

draw_proposals <- function(n = 1, mu, sigma){
  matrix(rnorm(n = length(mu) * n, mean = mu, sd = sigma), nrow = length(mu))
}
