library(devtools)
devtools::load_all()
mu <- 1:2
sigma <- 1:2
# sigma <- rep(0.001, 2)
logposterior <- lposterior_2
xprop <- draw_proposals(n = 3, mu, sigma)
eval_proposals(xprop, mu, sigma)
eval_logposterior(xprop, logposterior = logposterior)
gen_chain_apis(T = 1, M = 2, mu, sigma, logposterior = logposterior)
set.seed(1)
apis_chain <- gen_chain_apis(T = 100, M = 3, mu, sigma, logposterior = logposterior)
pmc_chain <- gen_chain_pmc(T = 100, M = 3, mu, sigma, logposterior = logposterior)
mcmc_chain <- gen_chain_mcmc(T = 100, M = 3, mu = c(10, -10), sigma, logposterior = logposterior, sigma_mh =rep(5, 2))
apis_pchain <- indep_chains_apis(d = 2, N = 20, T = 100, M = 3, mu, sigma, logposterior = logposterior)
pmc_pchain <- indep_chains_pmc(d = 2, N = 20, T = 100, M = 3, mu, sigma, logposterior = logposterior)
mcmc_pchain <- indep_chains_mcmc(d = 2, N =20, T = 100, M = 3, mu, sigma, logposterior = logposterior, sigma_mh = rep(5, 2))
denom_bybox <- with(mcmc_pchain, compute_denomtable_bybox(x, mu, sigma))
denom_byrow <- with(mcmc_pchain, compute_denomtable_byrow(x, mu, sigma))
denom_bytable <- with(mcmc_pchain, compute_denomtable_bytable(x, mu, sigma))
denom_byrow <- with(apis_pchain, compute_denomtable_byrow(x, mu, sigma))
denom_byrow <- with(pmc_pchain, compute_denomtable_byrow(x, mu, sigma))
xmat <- xarr4d_tomatrix(mcmc_pchain$x)
xmat <- xarr4d_tomatrix(pmc_pchain$x)
xmat <- xarr4d_tomatrix(apis_pchain$x)

N <- 100
T <- 10000
d <- 2
lais_chain <- lais(d = d, N = N, T = T, M = 2, mu = rnorm(N * T * d, sd = 10), sigma = rep(sqrt(3), d), lposterior_1, sigma_mh = rep(1, d), compute_denom = compute_denomtable_byrow)
compute_expectation(lais_chain$x, lais_chain$w)
compute_expectation(lais_chain$x, lais_chain$pi)

pmc_chain <- indep_pmc(d = d, N = N, T = T, M = 2, mu = rnorm(N * T * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_1, compute_denom = compute_denomtable_byrow, reuse_weights = FALSE)
# pmc_chain <- pmc(d = d, N = N, T = T, M = 2, mu = rnorm(N * T * d, sd = 10) , sigma = rep(sqrt(13), d), lposterior_1, compute_denom = compute_denomtable_byrow, reuse_weights = FALSE)
compute_expectation(pmc_chain$x, pmc_chain$w)
compute_expectation(pmc_chain$x, pmc_chain$pi)

apis_chain <- indep_apis(d = d, N = N, T = T, M = 2, mu = rnorm(N * T * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_1, compute_denom = compute_denomtable_byrow)
# apis_chain <- apis(d = d, N = N, T = T, M = 2, mu = rnorm(N * T * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_2, compute_denom = compute_denomtable_byrow)
compute_expectation(apis_chain$x, apis_chain$w)
compute_expectation(apis_chain$x, apis_chain$pi)
# draw_proposals <- function(x, mu, sigma, M = 1){
#   vapply(seq.int(M), draw_proposal, numeric(M), USE.NAMES = FALSE)
# }
