library(memoise)
library(mvtnorm)
source("luca_posterior.R")
source("MH_layer_algo.R")

eval_proposal <- function(x, mu, sigma){
  exp(sum(dnorm(x = x, mean = mu, sd = sigma, log = TRUE)))
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

compute_apis_weights <- function(xmat, mu, sigma, pi){
  proposals <- eval_proposals(xmat, mu, sigma)
  weights <- exp(pi - log(proposals))
  if(all(weights == 0)){
    weights <- rep(1, length(weights))
    warning("all weights are zero, replace by equiprobable weights")
  }
  weights / sum(weights)
}

compute_apis_weights <- function(pi, denom){
  weights <- exp(pi - denom)
  if(all(weights == 0)){
    weights <- rep(1, length(weights))
    warning("all weights are zero, replace by equiprobable weights")
  }
  weights / sum(weights)
}


compute_nextmu_apis <- function(xmat, weights){
  #apply(sweep(xmat, 2, weights, "*"), 1, sum)
  apply(xmat, 1, weighted.mean, w = weights)
}

compute_nextmu_pmc <- function(xmat, weights){
  if(any(is.na(weights))) browser()
  xmat[, rmultinom(1, 1, weights)]
}

compute_nextstep_apis <- function(M = 2, mu, sigma, logposterior){
  xprop <- draw_proposals(n = M, mu, sigma)
  weights <- compute_apis_weights(xprop, mu, sigma, logposterior)
  nextmu <- compute_nextmu_apis(xprop, weights)
  list("mu" = nextmu, "x" = xprop)
}

create_chain_stept <- function(compute_denom, compute_nextmu){
  chain_stept <- function(d, M = 2, N = 20, mu, sigma, logposterior){
    mu_next <- array(dim = c(1, d, N))
    mu_cur <- array(mu, dim = c(1, d, N))
    pi_cur <- array(dim = c(1, M, N))
    denom_cur <- array(dim = c(1, M, N))
    x_cur <- aperm(array(apply(mu_cur, 2, draw_proposals, n = M, sigma = sigma), c(1, N, M, d)), c(1, 4, 3, 2))
    pi_cur[1,,] <- apply(x_cur, 3:4, logposterior)
    denom_cur[1,,] <- compute_denom(x_arr = x_cur, mu_arr = mu_cur, sigma = sigma)
    weights <- apply(exp(pi_cur - denom_cur), c(3), function(w){
      if(all(w == 0)) w <- rep(1, length(w))
      return(w / sum(w))
    })
    for(n in seq.int(N)){
      mu_next[1,,n] <- compute_nextmu(x_cur[1,,, n], weights[, n])
    }
    list(mu = mu_next, x = x_cur, pi = pi_cur, denom = denom_cur)
  }
}

apis_stept <- create_chain_stept(compute_denomtable_bybox, compute_nextmu_apis)
apis_stept(d = 2, M = 2, N = 20, mu = 1:40 , sigma = rep(0.001, 2), logposterior = lposterior_2) 
pmc_stept <- create_chain_stept(compute_denomtable_byrow, compute_nextmu_pmc)
pmc_stept(d = 2, M = 2, N = 20, mu = 1:40, sigma = rep(0.001, 2), logposterior = lposterior_2) 

create_full_chain <- function(chain_stept){    
  full_chain <- function(d, T = 100, M = 2, N = 20, mu, sigma, logposterior, compute_denom, compute_nextmu){
    mu_res <- array(dim = c(T, d, N)) 
    x_res <- array(dim = c(T, d, M, N)) 
    pi_res <- array(dim = c(T, M, N)) 
    denom_res <- array(dim = c(T, M, N)) 
    mu_res[1,,] <- mu 
    for(t in seq.int(T)){
      t_step <- chain_stept(d = d, M = M, N = N, mu = mu_res[t,,], sigma = sigma, logposterior = logposterior)
      if(t < T) mu_res[t + 1,,] <- t_step$mu
      x_res[t,,,] <- t_step$x
      pi_res[t,,] <- t_step$pi
      denom_res[t,,] <- t_step$denom
    }
    list(mu = mu_res, x = x_res, pi = pi_res, denom = denom_res)
  }
}

apis_fchain <- create_full_chain(apis_stept)
apis_fchain(d = 2, T = 10, M = 2, N = 20, mu = rnorm(2*20), sigma = rep(sqrt(13), 2), logposterior = lposterior_2) 
pmc_fchain <- create_full_chain(pmc_stept)
pmc_fchain(d = 2, T = 100, M = 10, N = 200, mu = rnorm(2*200), sigma = rep(sqrt(13), 2), logposterior = lposterior_2) 

create_gen_chain <- function(compute_nextmu){
  gen_chain <- function(T = 100, M = 2,  mu, sigma, logposterior){
    mu_res <- matrix(nrow = T, ncol = length(mu))
    x_res <- array(dim = c(T, length(mu), M))
    pi_res <- array(dim = c(T, M))
    mu_res[1, ] <- mu
    x_res[1, ,] <- draw_proposals(n = M, mu, sigma)
    pi_res[1, ] <- eval_logposterior(x_res[1, ,], logposterior)
    for(t in seq.int(T)[-1]){
      x_cur <- x_res[t-1, ,]
      weights <- compute_apis_weights(x_cur, mu_res[t-1, ], sigma, pi_res[t-1, ])
      if(any(is.na(weights)))browser()
      mu_res[t, ] <- compute_nextmu(x_cur, weights)
      x_res[t, ,] <- draw_proposals(n = M, mu_res[t, ], sigma)
      pi_res[t, ] <- eval_logposterior(x_res[t, ,], logposterior)
    }
    list(mu = mu_res, x = x_res, pi = pi_res)
  }
}
gen_chain_apis <- create_gen_chain(compute_nextmu_apis)
gen_chain_pmc <- create_gen_chain(compute_nextmu_pmc)

gen_chain_mcmc <- function(T = 100, M = 2, mu, sigma, sigma_mh = sigma, logposterior){
  mu_res <- simple_MH(ndim = length(mu), log_posterior = logposterior, sigma = sigma_mh, init = mu, nsamples = T)
  x_res <- array(dim = c(T, length(mu), M))
  pi_res <- array(dim = c(T, M))
  for(t in seq.int(T)){
    x_res[t, ,] <- draw_proposals(n = M, mu_res[t, ], sigma)
    pi_res[t, ] <- eval_logposterior(x_res[t, ,], logposterior)
  }
  list(mu = mu_res, x = x_res, pi = pi_res)
}

make_parallel_chain <- function(gen_chain){
  parallel_chain <- function(d, N = 10, T = 100, M = 2, mu, sigma, logposterior, ...){
    init_mu <- matrix(mu, nrow = d, ncol = N)
    mu_res <- array(dim = c(T, d, N)) 
    x_res <- array(dim = c(T, d, M, N)) 
    pi_res <- array(dim = c(T, M, N)) 
    for(n in seq.int(N)){
      chain <- gen_chain(T = T, M = M, mu = init_mu[, n], sigma = sigma, logposterior = logposterior, ...)
      mu_res[,, n] <- chain$mu
      x_res[,,, n] <- chain$x
      pi_res[,, n] <- chain$pi
    }
    list(mu = mu_res, x = x_res, pi = pi_res)
  }
}
parallel_chain_mcmc <- make_parallel_chain(gen_chain_mcmc)
parallel_chain_pmc <- make_parallel_chain(gen_chain_pmc)
parallel_chain_apis <- make_parallel_chain(gen_chain_apis)

compute_phitable <- function(mu_arr, sigma){
  phi_res <- vector('list', length = dim(mu_arr)[1])
  for(i in seq_along(phi_res)){
#      phi_res[[i]] = apply(mu_arr[i,,, drop = FALSE], 2, function(mu){function(x) eval_proposal(x, mu =  mu, sigma = sigma)})
    for(j in seq.int(dim(mu_arr)[3])){
      phi_res[[i]][j] <-  list(function(x) eval_proposal(x, mu =  mu_arr[i,, j], sigma = sigma))
    }
  }
  phi_res
}

#compute_phitable <- function(mu_arr, sigma){
#  apply(mu_arr, c(1, 3), function(mu)list(function(x) eval_proposal(x, mu =  mu, sigma = sigma)))
#}

compute_pitable <- function(x_arr, logposterior){
  apply(x_arr, c(1, 3, 4), logposterior)
}

compute_denomtable_bybox <- function(x_arr, mu_arr, sigma){
  denom_res <- array(dim = dim(x_arr)[-2])
  for(i in seq.int(dim(denom_res)[1])){
    for(j in seq.int(dim(denom_res)[2])){
      for(k in seq.int(dim(denom_res)[3])){
        denom_res[i, j, k] <- eval_proposal(x_arr[i, ,j, k], mu_arr[i, ,k], sigma)
      }
    }
  }
  log(denom_res)
}

compute_denomtable_byrow <- function(x_arr, mu_arr, sigma){
  denom_res <- array(dim = dim(x_arr)[-2])
  for(i in seq.int(dim(denom_res)[1])){
    for(j in seq.int(dim(denom_res)[2])){
      for(k in seq.int(dim(denom_res)[3])){
        denom_res[i, j, k] <- sum(apply(mu_arr[i,,], 2, function(mu) eval_proposal(x_arr[i, ,j, k], mu, sigma)))
      }
    }
  }
  log(denom_res)
}

compute_denomtable_bytable <- function(x_arr, mu_arr, sigma){
  denom_res <- array(dim = dim(x_arr)[-2])
  for(i in seq.int(dim(denom_res)[1])){
    for(j in seq.int(dim(denom_res)[2])){
      for(k in seq.int(dim(denom_res)[3])){
        denom_res[i, j, k] <- sum(apply(apply(mu_arr, 2, identity), 1, function(mu) eval_proposal(x_arr[i, ,j, k], mu, sigma)))
      }
    }
  }
  log(denom_res)
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
  l_w <- wvect / sum(wvect)
  print(summary(l_w))
  apply(l_fx, 1, weighted.mean, w = l_w)
}

create_adaptive_is <- function(parallel_chain){
  adaptive_is <- function(d, N = 10, T = 100, M = 2, mu, sigma, logposterior, compute_denom = compute_denomtable_byrow, reuse_weights = FALSE, ...){
    pchain <- parallel_chain(d = d, N = N, T = T, M = M, mu = mu, sigma = sigma, logposterior = logposterior, ...)
    if(reuse_weights) denom_arr <- pchain$denom
    else denom_arr<- compute_denom(pchain$x, pchain$mu, sigma)
    weights <- compute_weights(pchain$pi, denom_arr)
    list(x = pchain$x, mu = pchain$mu, w = weights, pi = pchain$pi, denom = denom_arr, marglik = mean(weights))
  }
}  

lais <- create_adaptive_is(parallel_chain_mcmc)
apis <- create_adaptive_is(parallel_chain_apis)
apis <- create_adaptive_is(apis_fchain )
pmc <- create_adaptive_is(parallel_chain_pmc)
pmc <- create_adaptive_is(pmc_fchain)
        

mu <- 1:2
sigma <- 1:2
sigma <- rep(0.001, 2)
logposterior <- lposterior_2
xprop <- draw_proposals(n = 3, mu, sigma)
eval_proposals(xprop, mu, sigma)
eval_logposterior(xprop, logposterior = logposterior)
weights <- compute_apis_weights(xprop, mu, sigma, logposterior)
compute_nextmu_apis(xprop, weights)
compute_nextmu_pmc(xprop, weights)
compute_nextstep_apis(M = 2, mu, sigma, logposterior)
gen_chain_apis(T = 1, M = 2, mu, sigma, logposterior)
set.seed(1)
apis_chain <- gen_chain_apis(T = 100, M = 3, mu, sigma, logposterior)
pmc_chain <- gen_chain_pmc(T = 100, M = 3, mu, sigma, logposterior)
mcmc_chain <- gen_chain_mcmc(T = 10000, M = 3, mu = c(10, -10), sigma, lposterior_1, sigma_mh =rep(5, 2))
apis_pchain <- parallel_chain_apis(d = 2, N = 20, T = 100, M = 3, mu, sigma, logposterior)
pmc_pchain <- parallel_chain_pmc(d = 2, N = 20, T = 100, M = 3, mu, sigma, logposterior)
mcmc_pchain <- parallel_chain_mcmc(d = 2, N =20, T = 100, M = 3, mu, sigma, logposterior, sigma_mh = rep(5, 2))
phi_arr <- compute_phitable(mcmc_pchain$mu, sigma)
phi_arr <- compute_phitable(mcmc_pchain$mu[,, 1, drop = FALSE], sigma)
pi_arr <- compute_pitable(mcmc_pchain$x, logposterior)
pi_arr <- compute_pitable(mcmc_pchain$x[,,, 1, drop = FALSE], logposterior)
denom_bybox <- compute_denomtable_bybox(mcmc_pchain$x, phi_arr)
denom_byrow <- compute_denomtable_byrow(mcmc_pchain$x, phi_arr)
denom_byrow <- compute_denomtable_byrow(pmc_pchain$x, phi_arr)
denom_byrow <- compute_denomtable_byrow(apis_pchain$x, phi_arr)
denom_bytable <- compute_denomtable_bytable(mcmc_pchain$x, phi_arr)
w_byrow <- compute_weights(pi_arr, denom_byrow)
xmat <- xarr4d_tomatrix(mcmc_pchain$x)
xmat <- xarr4d_tomatrix(pmc_pchain$x)
xmat <- xarr4d_tomatrix(apis_pchain$x)
w_byrow_v <- as.vector(w_byrow)

lais_chain <- lais(d = 2, N = 20, T = 100, M = 100, mu, sigma = rep(sqrt(13), 2), logposterior, sigma_mh = rep(5, 2), compute_denom = compute_denomtable_byrow)
lais_chain <- lais(d = 2, N = 20, T = 100, M = 100, mu, sigma = rep(sqrt(13), 2), lposterior_1, sigma_mh = rep(5, 2), compute_denom = compute_denomtable_byrow)
lais_chain <- lais(d = 2, N = 20, T = 100, M = 100, mu, sigma = rep(sqrt(13), 2), lposterior_1, sigma_mh = rep(5, 2), compute_denom = compute_denomtable_bybox)
lais_chain <- lais(d = 2, N = 200, T = 100, M = 10, mu = rnorm(40, sd = 5), sigma = rep(sqrt(13), 2), lposterior_1, sigma_mh = rep(5, 2), compute_denom = compute_denomtable_byrow)
compute_expectation(lais_chain$x, lais_chain$w)
compute_expectation(lais_chain$x, lais_chain$pi)

pmc_chain <- pmc(d = 2, N = 20, T = 100, M = 100, mu, sigma = rep(sqrt(13), 2), logposterior, compute_denom = compute_denomtable_byrow)
pmc_chain <- pmc(d = 2, N = 20, T = 100, M = 100, mu, sigma = rep(sqrt(13), 2), logposterior, compute_denom = compute_denomtable_bybox)
pmc_chain <- pmc(d = 2, N = 200, T = 100, M = 10, mu = rnorm(400, sd = 5) , sigma = rep(sqrt(13), 2), lposterior_1, compute_denom = compute_denomtable_byrow, reuse_weights = TRUE)
compute_expectation(pmc_chain$x, pmc_chain$w)
compute_expectation(pmc_chain$x, pmc_chain$pi)

apis_chain <- apis(d = 2, N = 20, T = 100, M = 100, mu, sigma = rep(sqrt(13), 2), logposterior, compute_denom = compute_denomtable_byrow)
apis_chain <- apis(d = 2, N = 20, T = 100, M = 100, mu, sigma = rep(sqrt(13), 2), lposterior_1, compute_denom = compute_denomtable_byrow)
apis_chain <- apis(d = 2, N = 20, T = 100, M = 100, mu = rnorm(40, sd = 5) , sigma = rep(sqrt(13), 2), lposterior_1, compute_denom = compute_denomtable_bybox)
apis_chain <- apis(d = 2, N = 100, T = 100, M = 20, mu = rnorm(200, sd = 5) , sigma = rep(sqrt(13), 2), lposterior_1, compute_denom = compute_denomtable_byrow)
compute_expectation(apis_chain$x, apis_chain$w)
compute_expectation(apis_chain$x, apis_chain$pi)
# draw_proposals <- function(x, mu, sigma, M = 1){
#   vapply(seq.int(M), draw_proposal, numeric(M), USE.NAMES = FALSE)
# }