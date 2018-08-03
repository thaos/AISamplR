compute_nextmu_apis <- function(xmat, weights){
  #apply(sweep(xmat, 2, weights, "*"), 1, sum)
  if(any(is.na(weights))) stop("some weights are NA")
  apply(xmat, 1, weighted.mean, w = weights)
}

compute_nextmu_pmc <- function(xmat, weights){
  if(any(is.na(weights))) stop("some weights are NA")
  xmat[, as.logical(rmultinom(1, 1, weights))]
}

create_gen_mu_chain <- function(compute_nextmu){
  gen_mu_chain <- function(logposterior, mu, sigma, T = 100, M = 2){
    D  <- length(mu)
    mu_res <- matrix(nrow =  D, ncol = T)
    x_res <- array(dim = c(D, T, M))
    pi_res <- array(dim = c(T, M))
    denom_res <- array(dim = c(T, M))
    mu_res[, 1] <- mu
    x_res[, 1,] <- draw_proposals(n = M, mu, sigma)
    pi_res[1, ] <- eval_logposterior(x_res[, 1,], logposterior)
    denom_res[1, ] <- log(eval_proposals(x_res[, 1,], mu, sigma))
    for(t in seq.int(T)[-1]){
      weights <- normalize_weights(compute_weights(pi_res[t-1, ], denom_res[t-1, ]))
      mu_res[, t] <- compute_nextmu(x_res[, t-1,], weights)
      x_res[,t ,] <- draw_proposals(n = M, mu_res[, t], sigma)
      pi_res[t, ] <- eval_logposterior(x_res[, t,], logposterior)
      denom_res[t, ] <- log(eval_proposals(x_res[, t,], mu_res[, t], sigma))
    }
    return(mu_res)
  }
}

gen_mu_chain_pmc <- function(logposterior, mu, sigma, T = 100, M = 2){
  UseMethod("gen_mu_chain_pmc", logposterior)
}
gen_mu_chain_apis <- function(logposterior, mu, sigma, T = 100, M = 2){
  UseMethod("gen_mu_chain_apis", logposterior)
}
gen_mu_chain_apis.function <- create_gen_mu_chain(compute_nextmu_apis)
gen_mu_chain_pmc.function <- create_gen_mu_chain(compute_nextmu_pmc)

gen_mu_chain_pmc.externalptr <- function(logposterior, mu, sigma, T = 100, M = 2){
  gen_mu_chain_pmc_rcpp(
    lp = logposterior,
    mu = mu,
    sigma = sigma,
    T = T,
    M = M
  )
}
gen_mu_chain_apis.externalptr <- function(logposterior, mu, sigma, T = 100, M = 2){
  gen_mu_chain_apis_rcpp(
    lp = logposterior,
    mu = mu,
    sigma = sigma,
    T = T,
    M = M
  )
}

