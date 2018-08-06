// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "AISamplR_types.h"

#include <algorithm>    // std::min

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// from http://gallery.rcpp.org/articles/dmvnorm_arma/
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
    // Rcout << "proposition : "  << proposition << "; loglik : " << cur_lposterior <<std::endl ;
    rho = exp(cur_lposterior - prev_lposterior);
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

// [[Rcpp::export]]
NumericVector compute_denom_table_bybox_rcpp(NumericVector x, NumericVector mu, NumericVector sigma2, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector denom_res(T * N * M);
  denom_res.attr("dim") = IntegerVector::create(T, N, M);
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
    // Rcout << "denom : "  << sum(dmvnorm_arma(xmat, mumat, sigma_mat, true)) << std::endl ;
    denom_res[i] =  sum(dmvnorm_arma(xmat, mumat, sigma_mat, false));
  }  
  denom_res = log(denom_res);
  return(denom_res);
}

// [[Rcpp::export]]
NumericVector compute_denom_table_byrow_rcpp(NumericVector x, NumericVector mu, NumericVector sigma2, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector denom_res(T * N * M);
  denom_res.attr("dim") = IntegerVector::create(T, N, M);
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
  double denom = 0;
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
        denom = 0;
        for(std::size_t nprime = 0 ; nprime < N ; ++nprime){
          j = t + nprime * T;
          idx_mu = j *  D + (Rcpp::seq_len(D) - 1);
          for(std::size_t d = 0 ; d < D; ++d ){
            musub = mu[idx_mu];
            mumat(d) = musub(d);
          }
          denom += mean(dmvnorm_arma(xmat, mumat, sigma_mat, false));
        }
        denom_res[i] =  denom;
        p.increment(); // update progress
      }
    }
  }
  denom_res = log(denom_res);
  return(denom_res);
}

// [[Rcpp::export]]
NumericVector compute_denom_table_bytable_rcpp(NumericVector x, NumericVector mu, NumericVector sigma2, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector denom_res(T * N * M);
  denom_res.attr("dim") = IntegerVector::create(T, N, M);
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
  double denom = 0;
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
        denom = 0;
        for(std::size_t nprime = 0 ; nprime < N ; ++nprime){
          for(std::size_t tprime = 0 ; tprime < T ; ++tprime ){
            j = tprime + nprime * T;
            idx_mu = j *  D + (Rcpp::seq_len(D) - 1);
            for(std::size_t d = 0 ; d < D; ++d ){
              musub = mu[idx_mu];
              mumat(d) = musub(d);
            }
            denom += mean(dmvnorm_arma(xmat, mumat, sigma_mat, false));
          }
        }
        denom_res[i] =  denom;
        p.increment(); // update progress
      }
    }
  }
  denom_res = log(denom_res);
  return(denom_res);
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
  NumericVector denom_curr(M);
  NumericVector weight_curr(M);
  for(std::size_t t = 0 ; t < T - 1; ++t ){
    idx_mu = t *  D  + (Rcpp::seq_len(D) - 1);
    mu_curr = mu_chain[idx_mu];
    xs_curr = gen_xs_rcpp(mu_curr, sigma2, D, 1, 1, M);
    loglik_curr = compute_loglik_table_rcpp(lp, xs_curr, D, 1, 1, M);
    denom_curr = compute_denom_table_bybox_rcpp(xs_curr, mu_curr, sigma2, D, 1, 1, M);
    weight_curr = exp(loglik_curr - denom_curr);
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
  NumericVector denom_curr(M);
  NumericVector weight_curr(M);
  arma::ivec rn(M);
  for(std::size_t t = 0 ; t < T - 1; ++t ){
    idx_mu = t *  D  + (Rcpp::seq_len(D) - 1);
    mu_curr = mu_chain[idx_mu];
    xs_curr = gen_xs_rcpp(mu_curr, sigma2, D, 1, 1, M);
    loglik_curr = compute_loglik_table_rcpp(lp, xs_curr, D, 1, 1, M);
    denom_curr = compute_denom_table_bybox_rcpp(xs_curr, mu_curr, sigma2, D, 1, 1, M);
    weight_curr = exp(loglik_curr - denom_curr);
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

