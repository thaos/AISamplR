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

D <- 2

pmc_lp6_rcpp <- 
  pmc(lp6$pointer,
     mu = matrix(1:4, nrow = 2, ncol = 100),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = 100, T = 1000, M = 5)

with(pmc_lp6_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp6_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp6_rcpp, compute_expectation(x, weight))

pmc_lp6_r <- 
  pmc(lposterior_6,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = 100, T = 1000, M = 5)

with(pmc_lp6_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp6_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp6_r, compute_expectation(x, weight))

apis_lp6_rcpp <- 
  apis(lp6$pointer,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = 100, T = 1000, M = 5)

with(apis_lp6_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp6_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp6_rcpp, compute_expectation(x, weight))

apis_lp6_r <- 
  apis(lposterior_6,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = 100, T = 1000, M = 5)

with(apis_lp6_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp6_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp6_r, compute_expectation(x, weight))

lais_lp6_rcpp <- 
  lais(lp6$pointer,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = 100, T = 1000, M = 5)

with(lais_lp6_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp6_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp6_rcpp, compute_expectation(x, weight))

lais_lp6_r <- 
  lais(lposterior_6,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = 100, T = 1000, M = 5)

with(lais_lp6_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp6_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp6_r, compute_expectation(x, weight))
