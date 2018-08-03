create_gen_mu_pchains <- function(gen_mu_chain){
  gen_mu_chains <- function(logposterior,  mu, sigma2, T = 100, N = 2, ...){
    D <- nrow(mu)
    if(D != length(sigma2)) stop("different dimension for sigma2 and mu");
    if(N != ncol(mu)) stop("different dimension for N and nrow(mu)");
    mus_res <- array(NA, dim = c(D, T, N))
    for(n in seq.int(N)){
      mus_res[,, n] <- gen_mu_chain(
        logposterior = logposterior,
        mu = mu[, n],
        sigma2 = sigma2,
        T = T,
        ...
      )
    }
    return(mus_res)
  }
}

gen_mu_chains_mcmc <- create_gen_mu_pchains(simple_MH)
gen_mu_chains_pmc <- create_gen_mu_pchains(gen_mu_chain_pmc)
gen_mu_chains_apis <- create_gen_mu_pchains(gen_mu_chain_apis)

