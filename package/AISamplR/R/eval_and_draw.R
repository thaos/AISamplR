eval_proposal <- function(x, mu, sigma2){
  ans <- exp(sum(dnorm(x = x, mean = mu, sd = sqrt(sigma2), log = TRUE)))# + .Machine$double.xmin
  #if(any(is.na(ans) | ans == 0)) browser()
  ans
}

eval_proposals <- function(xmat, mu, sigma2){
  xmat <- matrix(xmat, nrow = nrow(xmat))
  apply(xmat, 2, eval_proposal, mu = mu, sigma2 = sigma2)
}

eval_logposterior <- function(xmat, logposterior){
  # dim(xmat) = c(D, 1, N)
  xmat <- matrix(xmat, nrow = nrow(xmat))
  apply(xmat, 2, logposterior)
}

draw_proposals <- function(n = 1, mu, sigma2){
  matrix(rnorm(n = length(mu) * n, mean = mu, sd = sqrt(sigma2)), nrow = length(mu))
}
