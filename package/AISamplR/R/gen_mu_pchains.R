create_gen_mu_pchains <- function(gen_mu_chain){
  gen_mu_chains <- function(logposterior,  mu, sigma, T = 100, N = 2, ...){
    D <- nrow(mu)
    if(D != length(sigma)) stop("different dimension for sigma and mu");
    if(N != ncol(mu)) stop("different dimension for N and nrow(mu)");
    mus_res <- array(NA, dim = c(D, T, N))
    for(n in seq.int(N)){
      mus_res[,, n] <- gen_mu_chain(
        logposterior = logposterior,
        mu = mu[n, ],
        sigma = sigma,
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

# logtarget <-  -1/32 * (4 - 10 * x1 - x2^2)^2 - x1^2/50 - x2^2/50
# body_lp6 <- '
# // [[Rcpp::export]]
#   double lposterior(NumericVector x){
#   double x1 = x[0];
#   double x2 = x[1];
#   double logtarget = -1.0/32 * pow(4 - 10 * x1 - pow(x2, 2), 2) - pow(x1, 2)/50 - pow(x2, 2)/50;
#   return logtarget;
# }
# '  
# lp6 <- make_lposterior_rcpp(body = body_lp6)
# 
# system.time(
#   gen_mu_chains_mcmc(
#     lp6$pointer,
#     mu = matrix(1:4, ncol = 2, nrow = 100),
#     sigma = c(1,1),
#     T = 1000, N = 100)
# )
# system.time(
#   gen_mu_chains_mcmc(
#     lp6$fun,
#     mu = matrix(1:4, ncol = 2, nrow = 100),
#     sigma = c(1,1),
#     T = 1000, N = 100
#   )
# )
# system.time(
#   gen_mu_chains_mcmc(
#     lposterior_6,
#     mu = matrix(1:4, ncol = 2, nrow = 100),
#     sigma = c(1,1),
#     T = 1000, N = 100
#   )
# )
# simple_MH(logposterior = lposterior_6,
#           mu = 1:2,
#           sigma = 1:2,
#           T = 100)
# simple_MH(logposterior = lp6$pointer,
#           mu = 1:2,
#           sigma = 1:2,
#           T = 100)
# 
gen_mu_chains_pmc(
                   lp6$pointer,
                   mu = matrix(1:4, nrow = 2, ncol = 2),
                   sigma = c(1,1),
                   T = 100, N = 2, M = 2 
                   )
