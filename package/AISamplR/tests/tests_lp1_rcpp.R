# mu_true = rep(1.6, 1.4);
body_lp1 <- '
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
  double muarr_1[2] = {-10, -10};
  double muarr_2[2] = {0, 16};
  double muarr_3[2] = {13, 8};
  double muarr_4[2] = {-9, 7};
  double muarr_5[2] = {14, -14};
  arma::rowvec mu1(&muarr_1[0], 2);
  arma::rowvec mu2(&muarr_2[0], 2);
  arma::rowvec mu3(&muarr_3[0], 2);
  arma::rowvec mu4(&muarr_4[0], 2);
  arma::rowvec mu5(&muarr_5[0], 2);
  double sigarr_1[4] = {2, 0.6, 0.6, 1};
  double sigarr_2[4] = {2, -0.4, -0.4, 2};
  double sigarr_3[4] = {2, 0.8, 0.8, 2};
  double sigarr_4[4] = {3, 0, 0, 5};
  double sigarr_5[4] = {2, -0.1, -0.1, 2};
  arma::mat sigma1(&sigarr_1[0], 2, 2);
  arma::mat sigma2(&sigarr_2[0], 2, 2);
  arma::mat sigma3(&sigarr_3[0], 2, 2);
  arma::mat sigma4(&sigarr_4[0], 2, 2);
  arma::mat sigma5(&sigarr_5[0], 2, 2);
  arma::mat xmat(1, x.length());
  for(std::size_t i = 0 ; i < x.length() ; ++i ){
    xmat(0, i) = x(i); 
  }
  double f_1 = sum(dmvnorm_arma(xmat, mu1, sigma1, false));
  double f_2 = sum(dmvnorm_arma(xmat, mu2, sigma2, false));
  double f_3 = sum(dmvnorm_arma(xmat, mu3, sigma3, false));
  double f_4 = sum(dmvnorm_arma(xmat, mu4, sigma4, false));
  double f_5 = sum(dmvnorm_arma(xmat, mu5, sigma5, false));
  // Rcout << "f_1 : "  << f_1 << std::endl ;
  // Rcout << "f_2 : "  << f_2 << std::endl ;
  // Rcout << "f_3 : "  << f_3 << std::endl ;
  // Rcout << "f_4 : "  << f_4 << std::endl ;
  // Rcout << "f_5 : "  << f_5 << std::endl ;
  double f = f_1 + f_2 + f_3 + f_4 + f_5;
  return log(f) - log(5);
}
'  
lp1 <- make_lposterior_rcpp(body = body_lp1)

D <- 2
T <- 100
N <- 20
M <- 5

pmc_lp1_rcpp <- 
  pmc(lp1$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(pmc_lp1_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp1_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp1_rcpp, compute_expectation(x, weight))

pmc_lp1_r <- 
  pmc(lposterior_1,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(pmc_lp1_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(pmc_lp1_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(pmc_lp1_r, compute_expectation(x, weight))

apis_lp1_rcpp <- 
  apis(lp1$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(apis_lp1_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp1_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp1_rcpp, compute_expectation(x, weight))

apis_lp1_r <- 
  apis(lposterior_1,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(apis_lp1_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(apis_lp1_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(apis_lp1_r, compute_expectation(x, weight))

lais_lp1_rcpp <- 
  lais(lp1$pointer,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(lais_lp1_rcpp, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp1_rcpp, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp1_rcpp, compute_expectation(x, weight))

lais_lp1_r <- 
  lais(lposterior_1,
     mu = matrix(1:4, nrow = D, ncol = N),
     sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
     compute_logdenom = compute_logdenom_byrow,
     N = N, T = T, M = M)

with(lais_lp1_r, rgl_plot(x[1,,,], x[2,,,], exp(loglik)))
with(lais_lp1_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
with(lais_lp1_r, compute_expectation(x, weight))
