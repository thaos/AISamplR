
# logtarget <-  -1/32 * (4 - 10 * x1 - x2^2)^2 - x1^2/50 - x2^2/50
body_lp6 <- '
// [[Rcpp::export]]
double lposterior(NumericVector x){
  double x1 = x[0];
  double x2 = x[1];
  double logtarget = -1.0/32 * pow(4 - 10 * x1 - pow(x2, 2), 2) - pow(x1, 2)/50 - pow(x2, 2)/50;
  return logtarget;
}
'  
lp6 <- make_lposterior_rcpp(body = body_lp6)

system.time(
            gen_mu_chains_mcmc(
                               lp6$pointer,
                               mu = matrix(1:4, nrow = 2, ncol = 100),
                               sigma = c(1,1),
                               T = 1000, N = 100)
            )
system.time(
            gen_mu_chains_mcmc(
                               lp6$fun,
                               mu = matrix(1:4, nrow = 2, ncol = 100),
                               sigma = c(1,1),
                               T = 1000, N = 100
                               )
            )
system.time(
            gen_mu_chains_mcmc(
                               lposterior_6,
                               mu = matrix(1:4, nrow = 2, ncol = 100),
                               sigma = c(1,1),
                               T = 1000, N = 100
                               )
            )
simple_MH(logposterior = lposterior_6,
          mu = 1:2,
          sigma = 1:2,
          T = 100)
simple_MH(logposterior = lp6$pointer,
          mu = 1:2,
          sigma = 1:2,
          T = 100)

gen_mu_chains_pmc(
                   lp6$pointer,
                   mu = matrix(1:4, nrow = 2, ncol = 2),
                   sigma = c(1,1),
                   T = 100, N = 2, M = 2 
                   )

gen_mu_chains_pmc(
                   lp6$fun,
                   mu = matrix(1:4, nrow = 2, ncol = 2),
                   sigma = c(1,1),
                   T = 100, N = 2, M = 2 
                   )

gen_mu_chains_apis(
                   lp6$pointer,
                   mu = matrix(1:4, nrow = 2, ncol = 2),
                   sigma = c(1,1),
                   T = 100, N = 2, M = 2 
                   )

gen_mu_chains_apis(
                   lp6$fun,
                   mu = matrix(1:4, nrow = 2, ncol = 2),
                   sigma = c(1,1),
                   T = 100, N = 2, M = 2 
                   )
