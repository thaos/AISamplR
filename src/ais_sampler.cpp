// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "AISamplR_types.h"

#include <algorithm>    // std::min
#include <math.h>       /* isinf, sqrt */


// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// from http://gallery.rcpp.org/articles/dmvnorm_arma/
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

// from http://gallery.rcpp.org/articles/simulate-multivariate-normal/
arma::mat rmvnorm_arma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
NumericMatrix simple_MH_rcpp(fp_logposterior_t lp, NumericVector mu, NumericVector sigma2, int T = 100){
  //Rcout << "simple_MH called \n";
  if(mu.length() != sigma2.length()) stop("different dimension for sigma2 and mu");
  NumericMatrix x(mu.length(), T);
  x(_, 0) = mu;
  arma::mat sigma_mat(mu.length(), mu.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < mu.length() ; ++i ){
    sigma_mat(i, i) = sigma2(i); 
  }
  // Rcout << "sigma_mat : "  << sigma_mat << std::endl ;
  NumericVector proposition;
  double prev_lposterior = (*lp)(x(_, 0));
  double cur_lposterior;
  double rho;
  double alpha;
  double u;
  for(std::size_t i = 1 ; i < T ; ++i ){
    // Rcout << "x : " << x(0, i - 1) << " ; "<< x(1, i - 1) << std::endl;
    proposition = as<NumericVector>(wrap(rmvnorm_arma(1, x(_, i - 1), sigma_mat)));
    cur_lposterior = (*lp)(proposition);
    if( std::isinf(cur_lposterior) && cur_lposterior < 0){
      rho = 0;
    } else {
    // Rcout << "proposition : "  << proposition << "; loglik : " << cur_lposterior <<std::endl ;
      rho = exp(cur_lposterior - prev_lposterior);
    }
    alpha = std::min(1.0, rho);
    u = arma::randu();
    if (u <= alpha){
      x(_, i) = proposition;
      prev_lposterior = cur_lposterior;
    } else{
      x(_, i) = x(_, i-1);
    }
  } 
  return x;
}

// [[Rcpp::export]]
NumericVector gen_xs_rcpp(NumericVector mu, NumericVector sigma2, int D, int T, int N, int M){
  if(D * T * N != mu.length()) stop("D * T * N != mu.length()");
  NumericVector x_res(D * T * N * M);
  x_res.attr("dim") = IntegerVector::create(D, T, N, M);
  int idx_x;
  int idx_mu;
  for(std::size_t m = 0 ; m < M ; ++m ){
    for(std::size_t n = 0 ; n < N ; ++n ){
      for(std::size_t t = 0 ; t < T ; ++t ){
        for(std::size_t d = 0 ; d < D ; ++d ){
          idx_x = d + t * D +  n * D * T + m * D * T * N;
          idx_mu = d + t * D +  n * D * T;
          x_res(idx_x) = arma::randn() * sqrt(sigma2(d)) + mu(idx_mu);
        }
      }
    }
  }
  return(x_res);
}

// [[Rcpp::export]]
NumericVector compute_loglik_table_rcpp(fp_logposterior_t lp, NumericVector x, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector loglik_res(T * N * M);
  loglik_res.attr("dim") = IntegerVector::create(T, N, M);
  IntegerVector idx;
  DoubleVector xsub(D); 
  for(std::size_t i = 0 ; i < T * N * M; ++i ){
    idx = i *  D  + (Rcpp::seq_len(D) - 1);
    // Rcout << "idx : "  << idx << std::endl ;
    xsub = x[idx];
    // Rcout << "xsub : "  << xsub << std::endl ;
    loglik_res[i] =  (*lp)(x[idx]);
  }  
  return(loglik_res);
}

//' @rdname compute_logdenom_bytable
//' @export
// [[Rcpp::export("compute_logdenom_bybox")]]
NumericVector compute_logdenom_bybox_rcpp(NumericVector x, NumericVector mu, NumericVector sigma2, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector logdenom_res(T * N * M);
  logdenom_res.attr("dim") = IntegerVector::create(T, N, M);
  IntegerVector idx_x;
  IntegerVector idx_mu;
  arma::mat xmat(1, D);
  arma::rowvec mumat(D);
  NumericVector xsub(D);
  NumericVector musub(D);
  arma::mat sigma_mat(sigma2.length(), sigma2.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < sigma2.length() ; ++i ){
    sigma_mat(i, i) = sigma2(i); 
  }
  for(std::size_t i = 0 ; i < T * N * M; ++i ){
    Rcpp::checkUserInterrupt() ;
    idx_x = i *  D  + (Rcpp::seq_len(D) - 1);
    // Rcout << "idx_x : "  << idx_x << std::endl ;
    idx_mu = i / M  * D + (Rcpp::seq_len(D) - 1);
    // Rcout << "idx_mu : "  << idx_mu << std::endl ;
    for(std::size_t d = 0 ; d < D; ++d ){
      xsub = x[idx_x];
      xmat(d) = xsub(d);
      musub = mu[idx_mu];
      mumat(d) = musub(d);
    }
    // Rcout << "i : "  << i << std::endl ;
    // Rcout << "xmat : "  << xmat << std::endl ;
    // Rcout << "mumat : "  << mumat << std::endl ;
    // Rcout << "sigma_mat : "  << sigma_mat << std::endl ;
    // Rcout << "logdenom : "  << sum(dmvnorm_arma(xmat, mumat, sigma_mat, true)) << std::endl ;
    logdenom_res[i] =  sum(dmvnorm_arma(xmat, mumat, sigma_mat, false));
  }  
  logdenom_res = log(logdenom_res);
  return(logdenom_res);
}

//' @rdname compute_logdenom_bytable
//' @export
// [[Rcpp::export("compute_logdenom_byrow")]]
NumericVector compute_logdenom_byrow_rcpp(NumericVector x, NumericVector mu, NumericVector sigma2, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector logdenom_res(T * N * M);
  logdenom_res.attr("dim") = IntegerVector::create(T, N, M);
  arma::mat sigma_mat(sigma2.length(), sigma2.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < sigma2.length() ; ++i ){
    sigma_mat(i, i) = sigma2(i); 
  }
  IntegerVector idx_x;
  IntegerVector idx_mu;
  arma::mat xmat(1, D);
  arma::rowvec mumat(D);
  NumericVector xsub(D);
  NumericVector musub(D);
  double logdenom = 0;
  int i;
  int j;
  Progress p(N * M * T, true);
  for(std::size_t m = 0 ; m < M ; ++m ){
    for(std::size_t n = 0 ; n < N ; ++n ){
      Rcpp::checkUserInterrupt() ;
      for(std::size_t t = 0 ; t < T ; ++t ){
        i = t + n * T + m * T * N;
        idx_x = i *  D  + (Rcpp::seq_len(D) - 1);
        xsub = x[idx_x];
        for(std::size_t d = 0 ; d < D; ++d ){
          xsub = x[idx_x];
          xmat(d) = xsub(d);
        }
        logdenom = 0;
        for(std::size_t nprime = 0 ; nprime < N ; ++nprime){
          j = t + nprime * T;
          idx_mu = j *  D + (Rcpp::seq_len(D) - 1);
          for(std::size_t d = 0 ; d < D; ++d ){
            musub = mu[idx_mu];
            mumat(d) = musub(d);
          }
          logdenom += mean(dmvnorm_arma(xmat, mumat, sigma_mat, false));
        }
        logdenom_res[i] =  logdenom;
        p.increment(); // update progress
      }
    }
  }
  logdenom_res = log(logdenom_res);
  return(logdenom_res);
}

//' Computes the logarithm of the denominator of importance sampling weights.
//' 
//' \code{compute_logdenom_bybox}, \code{compute_logdenom_byrow}
//' and \code{compute_logdenom_bytable}
//' offers three different ways to compute 
//' the logarithm of the denominator of the importance sampling weights.
//' 
//' The weighting step can be performed in different ways depending on 
//' how we considere the samples drawn at the sampling step.
//' We can think that one sample was drawn from its own corresponding proposal distribution 
//' or that one sample was drawn from a mixture of proposal distributions. 
//' In this package, three methods to compute the denominator of the weights are available:
//' \itemize{
//'   \item \code{compute_logdenom_bybox} 
//'   where we considered that a sample x_\{t,n,m\} 
//'   is drawn from an unique 
//'   proposal distribution with location parameter mu_\{t,n\}.
//'   \item \code{compute_logdenom_byrow} 
//'   where we considered that a sample  x_\{t,n,m\} 
//'   is drawn from an equiprobable mixture
//'   of all proposal distributions at time t.
//'   \item \code{compute_logdenom_bytable} 
//'   where we considered that a sample x_\{t,n,m\} 
//'   is drawn from an equiprobable mixture
//'   of all available proposal distributions, for all times t = 1,..., T and 
//'   chains of proposal n = 1,..., N.
//' }
//' 
//' @param x A numerical array of dimension D x T x N x M. 
//' It is the result of the sampling step in the adapative importance sampling scheme.
//' @param mu A numerical array of dimension D x T x N.
//' It is the result of the adaptaion step in the adapative importance sampling scheme.
//' @param sigma2  A numerical matrix of length D.
//' It provides the variance of each dimension of the gaussian proposal distribution
//' used in the sampling step.
//' @param D  An integer providing the dimension of the input space.
//' @param T  An integer providing the number iterations used in the adaptation step.
//' @param N  An integer providing the number of proposal distribution chains used.
//' @param M  An integer providing the number of samples drawn from each proposal.
//' 
//' @return An array of dimension T x N x M with the logarithm of the 
//' denominator needed to compute the importance sampling weights.
//' 
//' @examples
//' D <- 1
//' T <- 10
//' N <- 2
//' M <- 3
//' lpexp <- function(x){
//'   dexp(x, log = TRUE)
//' }
//' 
//' # Generates weighted samples for the exponential law
//' pmc_lpexp_r <- pmc(lpexp,
//'    mu = matrix(rnorm(D*N, mean = 1, sd = 1), nrow = D, ncol = N),
//'    sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
//'    compute_logdenom = compute_logdenom_byrow,
//'    N = N, T = T, M = M)
//' with(pmc_lpexp_r, plot(x  = x, y = weight))
//' # theorical value: ~ [1]
//' with(pmc_lpexp_r, compute_expectation(x, weight)) 
//' 
//' # Recompute the denominator
//' # with one gaussian proposal distribution associated to each sample.
//' logdenom_bybox <- 
//'   with(pmc_lpexp_r,
//'        compute_logdenom_bybox(x = x, mu = mu, sigma2 = 1,
//'                               D = D, T = T, N = N, M = M)
//'   )
//' # Plot new weights associated to each sample
//' with(pmc_lpexp_r, plot(x  = x, y = exp(loglik - logdenom_bybox)))
//' # Estimate expectation with the new weights
//' with(pmc_lpexp_r,
//'      compute_expectation(x  = x,
//'                          weight = exp(loglik - logdenom_bybox))
//' )
//' 
//' # Recompute the denominator
//' # with a proposal distribution which is a mixture of
//' # the N gaussian proposal distribution at time t.
//' logdenom_byrow <- 
//'  with(pmc_lpexp_r,
//'           compute_logdenom_byrow(x = x, mu = mu, sigma2 = 1,
//'                                  D = D, T = T, N = N, M = M)
//'  )
//' # Plot new weights associated to each sample
//' with(pmc_lpexp_r, plot(x  = x, y = exp(loglik - logdenom_byrow)))
//' # Estimate expectation with the new weights
//' with(pmc_lpexp_r,
//'      compute_expectation(x  = x,
//'                          weight = exp(loglik - logdenom_byrow))
//' )
//' 
//' # Recompute the denominator
//' # with a proposal distribution which is a mixture of
//' # the T x N gaussian proposal distribution available.
//' logdenom_bytable <- 
//'   with(pmc_lpexp_r,
//'        compute_logdenom_bytable(x = x, mu = mu, sigma2 = 1,
//'                                 D = D, T = T, N = N, M = M)
//'   )
//' # Plot new weights associated to each sample
//' with(pmc_lpexp_r, plot(x  = x, y = exp(loglik - logdenom_bytable)))
//' # Estimate expectation with the new weights
//' with(pmc_lpexp_r,
//'      compute_expectation(x  = x,
//'                          weight = exp(loglik - logdenom_bytable))
//' )
//' 
//' @export
// [[Rcpp::export("compute_logdenom_bytable")]]
NumericVector compute_logdenom_bytable_rcpp(NumericVector x, NumericVector mu, NumericVector sigma2, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector logdenom_res(T * N * M);
  logdenom_res.attr("dim") = IntegerVector::create(T, N, M);
  arma::mat sigma_mat(sigma2.length(), sigma2.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < sigma2.length() ; ++i ){
    sigma_mat(i, i) = sigma2(i); 
  }
  IntegerVector idx_x;
  IntegerVector idx_mu;
  arma::mat xmat(1, D);
  arma::rowvec mumat(D);
  NumericVector xsub(D);
  NumericVector musub(D);
  double logdenom = 0;
  int i;
  int j;
  Progress p(N * M * T, true);
  for(std::size_t m = 0 ; m < M ; ++m ){
    for(std::size_t n = 0 ; n < N ; ++n ){
      Rcpp::checkUserInterrupt() ;
      for(std::size_t t = 0 ; t < T ; ++t ){
        i = t + n * T + m * T * N;
        idx_x = i *  D  + (Rcpp::seq_len(D) - 1);
        xsub = x[idx_x];
        for(std::size_t d = 0 ; d < D; ++d ){
          xsub = x[idx_x];
          xmat(d) = xsub(d);
        }
        logdenom = 0;
        for(std::size_t nprime = 0 ; nprime < N ; ++nprime){
          for(std::size_t tprime = 0 ; tprime < T ; ++tprime ){
            j = tprime + nprime * T;
            idx_mu = j *  D + (Rcpp::seq_len(D) - 1);
            for(std::size_t d = 0 ; d < D; ++d ){
              musub = mu[idx_mu];
              mumat(d) = musub(d);
            }
            logdenom += mean(dmvnorm_arma(xmat, mumat, sigma_mat, false));
          }
        }
        logdenom_res[i] =  logdenom;
        p.increment(); // update progress
      }
    }
  }
  logdenom_res = log(logdenom_res);
  return(logdenom_res);
}

// [[Rcpp::export]]
NumericVector gen_mu_chain_apis_rcpp(fp_logposterior_t lp, NumericVector mu, NumericVector sigma2, int T, int M){
  int D = mu.length();
  // Rcout << "D: "  << D << std::endl ;
  NumericVector mu_chain(D * T);
  mu_chain.attr("dim") = IntegerVector::create(D, T);
  NumericVector xs_chain(D * T * M);
  xs_chain.attr("dim") = IntegerVector::create(D, T, M);
  arma::mat xt;
  IntegerVector idx_mu = (Rcpp::seq_len(D) - 1);
  int idx_x;
  mu_chain[idx_mu] = mu;
  NumericVector mu_curr(D);
  NumericVector xs_curr(D * M);
  NumericVector loglik_curr(M);
  NumericVector logdenom_curr(M);
  NumericVector weight_curr(M);
  for(std::size_t t = 0 ; t < T - 1; ++t ){
    idx_mu = t *  D  + (Rcpp::seq_len(D) - 1);
    mu_curr = mu_chain[idx_mu];
    xs_curr = gen_xs_rcpp(mu_curr, sigma2, D, 1, 1, M);
    loglik_curr = compute_loglik_table_rcpp(lp, xs_curr, D, 1, 1, M);
    logdenom_curr = compute_logdenom_bybox_rcpp(xs_curr, mu_curr, sigma2, D, 1, 1, M);
    weight_curr = exp(loglik_curr - logdenom_curr);
    bool zerotest = is_true(all(weight_curr == 0.0));
    if(zerotest) weight_curr = rep(1, weight_curr.length());
    weight_curr = weight_curr / sum(weight_curr);
    for(std::size_t d = 0 ; d < D ; ++d ){
      for(std::size_t m = 0 ; m < M ; ++m ){
        mu_chain(d + (t+1) * D) += weight_curr(m) * xs_curr(d + m * D);
        xs_chain(d + t * D + m * T * D) =  xs_curr(d + m * D);
      }
    }
  }
  return(mu_chain); 
}


// [[Rcpp::export]]
NumericVector gen_mu_chain_pmc_rcpp(fp_logposterior_t lp, NumericVector mu, NumericVector sigma2, int T, int M){
  int D = mu.length();
  // Rcout << "D: "  << D << std::endl ;
  NumericVector mu_chain(D * T);
  mu_chain.attr("dim") = IntegerVector::create(D, T);
  NumericVector xs_chain(D * T * M);
  xs_chain.attr("dim") = IntegerVector::create(D, T, M);
  arma::mat xt;
  IntegerVector idx_mu = (Rcpp::seq_len(D) - 1);
  mu_chain[idx_mu] = mu;
  NumericVector mu_curr(D);
  NumericVector xs_curr(D * M);
  NumericVector loglik_curr(M);
  NumericVector logdenom_curr(M);
  NumericVector weight_curr(M);
  arma::ivec rn(M);
  for(std::size_t t = 0 ; t < T - 1; ++t ){
    idx_mu = t *  D  + (Rcpp::seq_len(D) - 1);
    mu_curr = mu_chain[idx_mu];
    xs_curr = gen_xs_rcpp(mu_curr, sigma2, D, 1, 1, M);
    loglik_curr = compute_loglik_table_rcpp(lp, xs_curr, D, 1, 1, M);
    logdenom_curr = compute_logdenom_bybox_rcpp(xs_curr, mu_curr, sigma2, D, 1, 1, M);
    weight_curr = exp(loglik_curr - logdenom_curr);
    bool zerotest = is_true(all(weight_curr == 0.0));
    if(zerotest) weight_curr = rep(1, weight_curr.length());
    weight_curr = weight_curr / sum(weight_curr);
    for(std::size_t d = 0 ; d < D ; ++d ){
      rmultinom(1, weight_curr.begin() , M, rn.begin());
      for(std::size_t m = 0 ; m < M ; ++m ){
        if(rn(m) ==  1){
          mu_chain(d + (t+1) * D) = xs_curr(d + m * D);
        }
        xs_chain(d + t * D + m * T * D) =  xs_curr(d + m * D);
      }
    }
  }
  return(mu_chain); 
}

