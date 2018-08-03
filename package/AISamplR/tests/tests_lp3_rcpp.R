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
  double muarr_1[2] = {5, 0};
  double muarr_2[2] = {0, 16};
  arma::rowvec mu1(&muarr_1[0], 2);
  arma::rowvec mu2(&muarr_2[0], 2);
  double sigarr_1[4] = {2, 0.6, 0.6, 1};
  double sigarr_2[4] = {3, 0, 0, 3};
  arma::mat sigma1(&sigarr_1[0], d, d);
  arma::mat sigma2(&sigarr_2[0], d, d);
  arma::mat xmat(1, x.length());
  for(std::size_t i = 0 ; i < x.length() ; ++i ){
    xmat(0, i) = x(i); 
  }
  double f_1 = sum(dmvnorm_arma(xmat, mu1, sigma1, false));
  double f_2 = sum(dmvnorm_arma(xmat, mu2, sigma2, false));
  double f = f_1 + f_2;
  double ans = log(f) - log(2)
  return(ans);
}
'  
lp3 <- make_lposterior_rcpp(body = body_lp3)

D <- 2

pmc_lp3_rcpp <- 
  pmc(lp3$pointer,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(pmc_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp3_rcpp, compute_expectation(x, weight))

pmc_lp3_r <- 
  pmc(lposterior_3,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(pmc_lp3_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp3_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp3_r, compute_expectation(x, weight))

apis_lp3_rcpp <- 
  apis(lp3$pointer,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(apis_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp3_rcpp, compute_expectation(x, weight))

apis_lp3_r <- 
  apis(lposterior_3,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(apis_lp3_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp3_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp3_r, compute_expectation(x, weight))

lais_lp3_rcpp <- 
  lais(lp3$pointer,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(lais_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp3_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp3_rcpp, compute_expectation(x, weight))

lais_lp3_r <- 
  lais(lposterior_3,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = c(1, 1),sig2_prop = c(1, 1),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(lais_lp3_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp3_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp3_r, compute_expectation(x, weight))
