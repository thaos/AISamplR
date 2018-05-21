create_gen_chain <- function(compute_nextmu){
  gen_chain <- function(logposterior, T = 100, M = 2,  mu, sigma){
    mu_res <- matrix(nrow = T, ncol = length(mu))
    x_res <- array(dim = c(T, length(mu), M))
    pi_res <- array(dim = c(T, M))
    denom_res <- array(dim = c(T, M))
    mu_res[1, ] <- mu
    x_res[1, ,] <- draw_proposals(n = M, mu, sigma)
    pi_res[1, ] <- eval_logposterior(x_res[1, ,], logposterior)
    #denom_res[1, ] <- log(eval_proposals(x_res[1, ,], mu, sigma))
    for(t in seq.int(T)[-1]){
      weights <- normalize_weights(compute_weights(pi_res[t-1, ], denom_res[t-1, ]))
      mu_res[t, ] <- compute_nextmu(x_res[t-1, ,], weights)
      x_res[t, ,] <- draw_proposals(n = M, mu_res[t, ], sigma)
      pi_res[t, ] <- eval_logposterior(x_res[t, ,], logposterior)
     # denom_res[t, ] <- log(eval_proposals(x_res[t, ,], mu_res[t, ], sigma))
    }
    list(mu = mu_res, x = x_res, pi = pi_res)#, denom = denom_res)
  }
}
gen_chain_apis <- create_gen_chain(compute_nextmu_apis)
gen_chain_pmc <- create_gen_chain(compute_nextmu_pmc)

gen_chain_mcmc <- function(T = 100, M = 2, mu, sigma, sigma_mh = sigma, logposterior){
  mu_res <- simple_MH(logposterior = logposterior, d = length(mu), T = T, mu = mu, sigma = sigma_mh)
  x_res <- array(dim = c(T, length(mu), M))
  pi_res <- array(dim = c(T, M))
  #denom_res <- array(dim = c(T, M))
  for(t in seq.int(T)){
    x_res[t, ,] <- draw_proposals(n = M, mu_res[t, ], sigma)
    #denom_res[t, ] <- eval_proposals(x_res[t, ,], mu_res[t, ], sigma)
    pi_res[t, ] <- eval_logposterior(x_res[t, ,], logposterior)
  }
  list(mu = mu_res, x = x_res, pi = pi_res)#, denom = denom_res)
}

make_parallel_chain <- function(gen_chain){
  parallel_chain <- function(logposterior, d, N = 10, T = 100, M = 2, mu, sigma, ...){
    init_mu <- matrix(mu, nrow = d, ncol = N)
    mu_res <- array(dim = c(T, d, N))
    x_res <- array(dim = c(T, d, M, N))
    pi_res <- array(dim = c(T, M, N))
    #denom_res <- array(dim = c(T, M, N))
    for(n in seq.int(N)){
      chain <- gen_chain(T = T, M = M, mu = init_mu[, n], sigma = sigma, logposterior = logposterior, ...)
      mu_res[,, n] <- chain$mu
      x_res[,,, n] <- chain$x
      pi_res[,, n] <- chain$pi
      #denom_res[,, n] <- chain$denom
    }
    list(mu = mu_res, x = x_res, pi = pi_res)#, denom = denom_res)
  }
}
indep_chains_mcmc <- make_parallel_chain(gen_chain_mcmc)
indep_chains_pmc <- make_parallel_chain(gen_chain_pmc)
indep_chains_apis <- make_parallel_chain(gen_chain_apis)
