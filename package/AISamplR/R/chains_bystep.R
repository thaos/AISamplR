create_chain_stept <- function(compute_denom, compute_nextmu){
  chain_stept <- function(logposterior, d, M = 2, N = 20, mu, sigma){
    mu_next <- array(dim = c(1, d, N))
    mu_cur <- array(mu, dim = c(1, d, N))
    pi_cur <- array(dim = c(1, M, N))
    denom_cur <- array(dim = c(1, M, N))
    x_cur <- aperm(array(apply(mu_cur, 2, draw_proposals, n = M, sigma = sigma), c(1, N, M, d)), c(1, 4, 3, 2))
    pi_cur[1,,] <- apply(x_cur, 3:4, logposterior)
    denom_cur[1,,] <- compute_denom(x_arr = x_cur, mu_arr = mu_cur, sigma = sigma)
    weights <- apply(compute_weights(pi_cur, denom_cur), c(3), normalize_weights)
    for(n in seq.int(N)){
      mu_next[1,,n] <- compute_nextmu(x_cur[1,,, n], weights[, n])
    }
    list(mu = mu_next, x = x_cur, pi = pi_cur, denom = denom_cur)
  }
}

stept_apis <- create_chain_stept(compute_denomtable_bybox, compute_nextmu_apis)
#apis_stept(d = 2, M = 2, N = 20, mu = 1:40 , sigma = rep(0.001, 2), logposterior = lposterior_2)
stept_pmc <- create_chain_stept(compute_denomtable_byrow, compute_nextmu_pmc)
#pmc_stept(d = 2, M = 2, N = 20, mu = 1:40, sigma = rep(0.001, 2), logposterior = lposterior_2)

create_full_chain <- function(chain_stept){
  full_chain <- function(logposterior, d, T = 100, M = 2, N = 20, mu, sigma, compute_denom, compute_nextmu){
    mu_res <- array(dim = c(T, d, N))
    x_res <- array(dim = c(T, d, M, N))
    pi_res <- array(dim = c(T, M, N))
    denom_res <- array(dim = c(T, M, N))
    mu_res[1,,] <- mu
    for(t in seq.int(T)){
      t_step <- chain_stept(logposterior = logposterior, d = d, M = M, N = N, mu = mu_res[t,,], sigma = sigma)
      if(t < T) mu_res[t + 1,,] <- t_step$mu
      x_res[t,,,] <- t_step$x
      pi_res[t,,] <- t_step$pi
      denom_res[t,,] <- t_step$denom
    }
    list(mu = mu_res, x = x_res, pi = pi_res, denom = denom_res)
  }
}

fchain_apis <- create_full_chain(stept_apis)
#apis_fchain(d = 2, T = 10, M = 2, N = 20, mu = rnorm(2*20), sigma = rep(sqrt(13), 2), logposterior = lposterior_2)
fchain_pmc <- create_full_chain(stept_pmc)
#pmc_fchain(d = 2, T = 100, M = 10, N = 200, mu = rnorm(2*200), sigma = rep(sqrt(13), 2), logposterior = lposterior_2)
