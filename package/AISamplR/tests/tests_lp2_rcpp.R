# mu_true = rep(5, 10);
body_lp2 <- '
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
  double mu_arr[2] = {0, 16};
  double sigma_arr[4] = {3, 0, 0, 3};
  arma::rowvec mu(&mu_arr[0], 2);
  arma::mat sigma(&sigma_arr[0], 2, 2);
  // Rcout << "x : " << x << std::endl ;
  arma::mat xmat(1, x.length());
  for(std::size_t i = 0 ; i < x.length() ; ++i ){
    xmat(0, i) = x(i); 
  }
  // Rcout << "xmat : "  << xmat << std::endl ;
  // Rcout << "Je bug pas encore" << std::endl ;
  double ans = sum(dmvnorm_arma(xmat, mu, sigma, true));
  // Rcout << "x.length : " << x.length() << std::endl ;
  // for(std::size_t i = 0 ; i < x.length() ; ++i ){   
  //   Rcout << "i : " << i << std::endl ;
  //   Rcout << "x : " << x(i) << std::endl ;
  //   ans +=  R::dnorm(x(i), 0.0, 1.0, true) ;
  // }
  // Rcout << "ans : " << ans << std::endl ;
  return(ans) ;
}
'  
lp2 <- make_lposterior_rcpp(body = body_lp2)

D <- 10

pmc_lp2_rcpp <- 
  pmc(lp2$pointer,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(pmc_lp2_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp2_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp2_rcpp, compute_expectation(x, weight))

pmc_lp2_r <- 
  pmc(lposterior_2,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(pmc_lp2_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp2_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp2_r, compute_expectation(x, weight))

apis_lp2_rcpp <- 
  apis(lp2$pointer,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(apis_lp2_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp2_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp2_rcpp, compute_expectation(x, weight))

apis_lp2_r <- 
  apis(lposterior_2,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(apis_lp2_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp2_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp2_r, compute_expectation(x, weight))

lais_lp2_rcpp <- 
  lais(lp2$pointer,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(lais_lp2_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp2_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp2_rcpp, compute_expectation(x, weight))

lais_lp2_r <- 
  lais(lposterior_2,
     mu = matrix(1:4, nrow = D, ncol = 100),
     sig2_adapt = c(1, 1),sig2_prop = c(1, 1),
     sig2_adapt = rep(1, D), sig2_prop = rep(1, D),
     compute_denom = compute_denom_table_byrow_rcpp,
     N = 100, T = 1000, M = 5)

with(lais_lp2_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp2_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp2_r, compute_expectation(x, weight))
