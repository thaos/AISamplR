Sys.setenv("R_TESTS" = "")
library(AISamplR)

lposterior_3 <- function(x){
  mu_1 <- c(5, 0)
  sigma_1 <- matrix(c(2, rep(0.6, 2), 1), ncol = 2, nrow = 2)
  mu_2 <- c(0, 16)
  sigma_2 <- matrix(c(3, rep(0, 2), 3), ncol = 2, nrow = 2)
  f_1  <- 1/2 * mvtnorm::dmvnorm(x, mean = mu_1 , sigma = sigma_1) 
  f_2  <- 1/2 * mvtnorm::dmvnorm(x, mean = mu_2 , sigma = sigma_2) 
  f <- f_1 + f_2 
  log(f)
}


# mu_true = rep(2.5, 8);
body_lp3 <- '
const double log2pi = std::log(2.0 * M_PI);
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}
arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false) { 
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;

  if (log) { 
    return(logretval);
  } else { 
    return(exp(logretval));
  }
}
// [[Rcpp::export]]
double lposterior(NumericVector x){
  int d = 2;
  int d2 = pow(d, 2);
  double muarr_1[2] = {5.0, 0.0};
  double muarr_2[2] = {0.0, 16.0};
  arma::rowvec mu1(&muarr_1[0], 2);
  arma::rowvec mu2(&muarr_2[0], 2);
  double sigarr_1[4] = {2.0, 0.6, 0.6, 1.0};
  double sigarr_2[4] = {3.0, 0.0, 0.0, 3.0};
  arma::mat sigma1(&sigarr_1[0], d, d);
  arma::mat sigma2(&sigarr_2[0], d, d);
  arma::mat xmat(1, x.length());
  for(std::size_t i = 0 ; i < x.length() ; ++i ){
    xmat(0, i) = x(i); 
  }
  double f_1 = sum(dmvnorm_arma(xmat, mu1, sigma1, false));
  double f_2 = sum(dmvnorm_arma(xmat, mu2, sigma2, false));
  double f = f_1 + f_2;
  double ans = log(f) - log(2);
  return(ans);
}
'  
lp3 <- make_lposterior_rcpp(body = body_lp3)

D <- 2
T <- 100
N <- 20
M <- 5

pmc_lp3_rcpp <- 
  pmc(lp3$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(pmc_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp3_rcpp, compute_expectation(x, weight))

pmc_lp3_r <- 
  pmc(lposterior_3,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(pmc_lp3_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp3_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp3_r, compute_expectation(x, weight))

apis_lp3_rcpp <- 
  apis(lp3$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(apis_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp3_rcpp, compute_expectation(x, weight))

apis_lp3_r <- 
  apis(lposterior_3,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(apis_lp3_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp3_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp3_r, compute_expectation(x, weight))

lais_lp3_rcpp <- 
  lais(lp3$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(lais_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp3_rcpp, compute_expectation(x, weight))

lais_lp3_r <- 
  lais(lposterior_3,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(lais_lp3_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp3_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp3_r, compute_expectation(x, weight))
