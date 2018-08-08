lposterior_0 <- function(x){
  dexp(x, log = TRUE)
}

body_lp0 <- '
// [[Rcpp::export]]
double lposterior(NumericVector x){
  double logtarget; 
  if(x[0] < 0){
    logtarget = log(0.0);
  } else{
    logtarget = -x[0];
  }
  return logtarget;
}
'  
lp0 <- make_lposterior_rcpp(body = body_lp0)

D <- 1
N <- 10
T <- 100
M <- 5

pmc_lp0_rcpp <- 
  pmc(lp0$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(pmc_lp0_rcpp, rgl_plot(x[1,,,], 0, exp(loglik)))
with(pmc_lp0_rcpp, rgl_plot(x[1,,,], 0, sqrt(weight)))
with(pmc_lp0_rcpp, compute_expectation(x, weight))

pmc_lp0_r <- 
  pmc(lposterior_0,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(pmc_lp0_r, rgl_plot(x[1,,,], 0, exp(loglik)))
with(pmc_lp0_r, rgl_plot(x[1,,,], 0, sqrt(weight)))
with(pmc_lp0_r, compute_expectation(x, weight))

apis_lp0_rcpp <- 
  apis(lp0$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(apis_lp0_rcpp, rgl_plot(x[1,,,], 0, exp(loglik)))
with(apis_lp0_rcpp, rgl_plot(x[1,,,], 0, sqrt(weight)))
with(apis_lp0_rcpp, compute_expectation(x, weight))

apis_lp0_r <- 
  apis(lposterior_0,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(apis_lp0_r, rgl_plot(x[1,,,], 0, exp(loglik)))
with(apis_lp0_r, rgl_plot(x[1,,,], 0, sqrt(weight)))
with(apis_lp0_r, compute_expectation(x, weight))

lais_lp0_rcpp <- 
  lais(lp0$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(lais_lp0_rcpp, rgl_plot(x[1,,,], 0, exp(loglik)))
with(lais_lp0_rcpp, rgl_plot(x[1,,,], 0, sqrt(weight)))
with(lais_lp0_rcpp, compute_expectation(x, weight))

lais_lp0_r <- 
  lais(lposterior_0,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(lais_lp0_r, rgl_plot(x[1,,,], 0, exp(loglik)))
with(lais_lp0_r, rgl_plot(x[1,,,], 0, sqrt(weight)))
with(lais_lp0_r, compute_expectation(x, weight))
