Sys.setenv("R_TESTS" = "")
library(AISamplR)

lposterior_5 <- function(x){
  mu <- rep(5, 10)
  sigma <- diag(4, 10)
  mvtnorm::dmvnorm(x, mean = mu , sigma = sigma, log = TRUE) 
}

# mu_true = rep(5, 10);
body_lp5 <- '
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
  int d = 10;
  int d2 = pow(d, 2);
  double mu_arr[d];
  for(std::size_t i = 0 ; i < d ; ++i ){
    mu_arr[i] = 5.0;
  }
  arma::rowvec mu(&mu_arr[0], d);
  double sigma_arr[d2];
  for(std::size_t i = 0 ; i < d2 ; ++i ){
    sigma_arr[i] = 0;
  }
  std::size_t j;
  for(std::size_t i = 0 ; i < d; ++i ){
    j = i * 10 + i;
    sigma_arr[j] = 1;
  }
  arma::mat sigma(&sigma_arr[0], d, d);
  // Rcout << "x : "  << x << std::endl; 
  // Rcout << "mu : "  << mu << std::endl; 
  // Rcout << "sigma :"  << std::endl; 
  // Rcout << sigma << std::endl; 
  arma::mat xmat(1, d);
  for(std::size_t i = 0 ; i < d ; ++i ){
    xmat(0, i) = x(i); 
  }
  double ans = sum(dmvnorm_arma(xmat, mu, sigma, true));
  return(ans);
}
'  
lp5 <- make_lposterior_rcpp(body = body_lp5)

D <- 10
T <- 100
N <- 20
M <- 5

pmc_lp5_rcpp <- 
  pmc(lp5$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(pmc_lp5_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp5_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp5_rcpp, compute_expectation(x, weight))

pmc_lp5_r <- 
  pmc(lposterior_5,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(pmc_lp5_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp5_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp5_r, compute_expectation(x, weight))

apis_lp5_rcpp <- 
  apis(lp5$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(apis_lp5_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp5_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp5_rcpp, compute_expectation(x, weight))

apis_lp5_r <- 
  apis(lposterior_5,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(apis_lp5_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp5_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp5_r, compute_expectation(x, weight))

lais_lp5_rcpp <- 
  lais(lp5$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(lais_lp5_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp5_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp5_rcpp, compute_expectation(x, weight))

lais_lp5_r <- 
  lais(lposterior_5,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(lais_lp5_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp5_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp5_r, compute_expectation(x, weight))
