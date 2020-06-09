compute_nextmu_apis <- function(xmat, weights){
  #apply(sweep(xmat, 2, weights, "*"), 1, sum)
  if(any(is.na(weights))) stop("some weights are NA")
  apply(xmat, 1, weighted.mean, w = weights)
}

compute_nextmu_pmc <- function(xmat, weights){
  if(any(is.na(weights))) stop("some weights are NA")
  xmat[, ,as.logical(rmultinom(1, 1, weights))]
}

create_gen_mu_chain <- function(compute_nextmu){
  gen_mu_chain <- function(logposterior, mu, sigma2, T = 100, M = 2){
    D  <- length(mu)
    mu_res <- matrix(nrow =  D, ncol = T)
    x_res <- array(dim = c(D, T, M))
    pi_res <- array(dim = c(T, M))
    logdenom_res <- array(dim = c(T, M))
    mu_res[, 1] <- mu
    x_res[, 1,] <- draw_proposals(n = M, mu, sigma2)
    pi_res[1, ] <- eval_logposterior(x_res[, 1,, drop = FALSE], logposterior)
    logdenom_res[1, ] <- log(eval_proposals(x_res[, 1,, drop = FALSE], mu, sigma2))
    for(t in seq.int(T)[-1]){
      weights <- normalize_weights(
        compute_weight_table(pi_res[t-1, ], logdenom_res[t-1, ])
      )
      mu_res[, t] <- compute_nextmu(x_res[, t-1,, drop = FALSE], weights)
      x_res[,t ,] <- draw_proposals(n = M, mu_res[, t, drop = FALSE], sigma2)
      pi_res[t, ] <- eval_logposterior(x_res[, t,, drop = FALSE], logposterior)
      logdenom_res[t, ] <- log(eval_proposals(x_res[, t,, drop = FALSE],
                                              mu_res[, t, drop = FALSE],
                                              sigma2))
    }
    return(mu_res)
  }
}

gen_mu_chain_pmc <- function(logposterior, mu, sigma2, T = 100, M = 2){
  UseMethod("gen_mu_chain_pmc", logposterior)
}
gen_mu_chain_apis <- function(logposterior, mu, sigma2, T = 100, M = 2){
  UseMethod("gen_mu_chain_apis", logposterior)
}
gen_mu_chain_apis.function <- create_gen_mu_chain(compute_nextmu_apis)
gen_mu_chain_pmc.function <- create_gen_mu_chain(compute_nextmu_pmc)

gen_mu_chain_pmc.externalptr <- function(logposterior, mu, sigma2, T = 100, M = 2){
  gen_mu_chain_pmc_rcpp(
    lp = logposterior,
    mu = mu,
    sigma2 = sigma2,
    T = T,
    M = M
  )
}
gen_mu_chain_apis.externalptr <- function(logposterior, mu, sigma2, T = 100, M = 2){
  gen_mu_chain_apis_rcpp(
    lp = logposterior,
    mu = mu,
    sigma2 = sigma2,
    T = T,
    M = M
  )
}

