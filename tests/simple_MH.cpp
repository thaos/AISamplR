
#include <algorithm>    // std::min
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
// [[Rccp::plugins(cpp11)]]
// [[Rcpp::export]]


NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

typedef double (*funcPtrG)(NumericVector x);
typedef XPtr<funcPtrG> fptr_t;

typedef double (*fp_logposterior)(NumericVector x);
typedef XPtr<fp_logposterior> fp_logposterior_t;

#include <RcppArmadillo.h>

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
double logposterior(NumericVector x){
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

// // [[Rcpp::export]]
// double GGfuncpp(NumericVector x)
// {
//   Rcout << "GGfuncpp called\n";
//   return sum(x);
// }
// 
// [[Rcpp::export]]
fp_logposterior_t make_logposteterior_cpp()
{
  return fp_logposterior_t(new fp_logposterior(logposterior));
}
// 
// 
// // [[Rcpp::export]]
// fptr_t make_GGfuncpp()
// {
//   return fptr_t(new funcPtrG(GGfuncpp));
// }
// 
// // [[Rcpp::export]]
// double testfun(fptr_t fun2, NumericVector x)
// {
//   return (*fun2)(x);
// }
// 
// 
// [[Rcpp::export]]
NumericMatrix simple_MH_rcpp(fp_logposterior_t lp, NumericVector mu, NumericVector sigma, int T = 100){
  Rcout << "simple_MH called \n";
  if(mu.length() != sigma.length()) stop("different dimension for sigma and mu");
  NumericMatrix x(T, mu.length());
  x(0, _) = mu;
  arma::mat sigma_mat(mu.length(), mu.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < mu.length() ; ++i ){
    sigma_mat(i, i) = sigma(i); 
  }
  Rcout << "sigma_mat : "  << sigma_mat << std::endl ;
  NumericVector proposition;
  double prev_lposterior = (*lp)(x(0, _));
  double cur_lposterior;
  double rho;
  double alpha;
  double u;
  for(std::size_t i = 1 ; i < T ; ++i ){
    // Rcout << "x : " << x(i - 1, 0) << " ; "<< x(i - 1 , 1) << std::endl;
    proposition = as<NumericVector>(wrap(rmvnorm_arma(1, x(i - 1, _), sigma_mat)));
    cur_lposterior = (*lp)(proposition);
    // Rcout << "proposition : "  << proposition << "; loglik : " << cur_lposterior <<std::endl ;
    rho = exp(cur_lposterior - prev_lposterior);
    alpha = std::min(1.0, rho);
    u = arma::randu();
    if (u <= alpha){
      x(i, _) = proposition;
      prev_lposterior = cur_lposterior;
    } else{
        x(i, _) = x(i-1, _);
    }
  } 
  return x;
}

// [[Rcpp::export]]
Rcpp::List gen_chain_mcmc_rcpp(fp_logposterior_t lp, NumericVector mu, NumericVector sigma, int T = 100, int M = 2){
  NumericMatrix mu_res = simple_MH_rcpp(lp, mu, sigma, T);
  NumericVector x_res(mu.length() * T * M);
  x_res.attr("dim") = IntegerVector::create(T, mu.length(), M);
  arma::mat sigma_mat(mu.length(), mu.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < mu.length() ; ++i ){
    sigma_mat(i, i) = sigma(i); 
  }
  int idx;
  arma::uvec pos(1);
  int D = mu.length();
  arma::mat xdrawn(D, M);
  for(std::size_t t = 0 ; t < T ; ++t ){
    xdrawn = rmvnorm_arma(M, mu_res(t, _), sigma_mat);
    for (std::size_t d = 0 ; d < D ; ++d ){
      for (std::size_t m = 0 ; m < M ; ++m ){
        idx = t + d * T + m * T * D; 
        Rcout << "d : "  << d <<std::endl ;
        Rcout << "m: "  << m <<std::endl ;
        Rcout << "idx : "  << idx <<std::endl ;
        Rcout << "xdrawn(d, m): "  << xdrawn(d, m) <<std::endl ;
        x_res(idx) = xdrawn(d, m);
      }
    }
  }
  return Rcpp::List::create(mu_res, x_res);
}

// [[Rcpp::export]]
NumericVector gen_mus_pmcmc_rcpp(fp_logposterior_t lp, NumericMatrix mu, NumericVector sigma, int T = 100, int N = 2){
  int D = mu.ncol();
  if(D != sigma.length()) stop("different dimension for sigma and mu");
  if(N != mu.nrow()) stop("different dimension for N and mu.nrow()");
  NumericVector mus_res(D * T * N);
  mus_res.attr("dim") = IntegerVector::create(D, T, N);
  int idx_start; 
  int idx_end;
  IntegerVector idx;
  for(std::size_t n = 0 ; n < N ; ++n ){
    idx_start = n * D * T;
    idx_end = (n + 1) * T * D - 1;
    idx = seq_len(idx_end - idx_start + 1 ) - 1 + idx_start;
    Rcout << "idx : "  << idx <<std::endl ;
    mus_res[idx] =  transpose(simple_MH_rcpp(lp, mu(n, _), sigma, T));
  }  
  return(mus_res);
}

// [[Rcpp::export]]
NumericVector gen_xs_rcpp(NumericVector mu, NumericVector sigma, int D, int T, int N, int M){
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
          x_res(idx_x) = arma::randn() * sqrt(sigma(d)) + mu(idx_mu);
        }
      }
    }
  }
  return(x_res);
}

// [[Rcpp::export]]
NumericVector create_loglik_table_rcpp(fp_logposterior_t lp, NumericVector x, int D, int T,  int N, int M){
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
NumericVector create_denom_table_bybox_rcpp(NumericVector x, NumericVector mu, NumericVector sigma, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector denom_res(T * N * M);
  denom_res.attr("dim") = IntegerVector::create(T, N, M);
  IntegerVector idx_x;
  IntegerVector idx_mu;
  arma::mat xmat(1, D);
  arma::rowvec mumat(D);
  NumericVector xsub(D);
  NumericVector musub(D);
  arma::mat sigma_mat(sigma.length(), sigma.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < sigma.length() ; ++i ){
    sigma_mat(i, i) = sigma(i); 
  }
  for(std::size_t i = 0 ; i < T * N * M; ++i ){
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
  return(log(denom_res));
}

// [[Rcpp::export]]
NumericVector create_denom_table_byrow_rcpp(NumericVector x, NumericVector mu, NumericVector sigma, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector denom_res(T * N * M);
  denom_res.attr("dim") = IntegerVector::create(T, N, M);
  arma::mat sigma_mat(sigma.length(), sigma.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < sigma.length() ; ++i ){
    sigma_mat(i, i) = sigma(i); 
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
  // or(std::size_t i = 0 ; i < T * N * M; ++i ){
  for(std::size_t m = 0 ; m < M ; ++m ){
    for(std::size_t n = 0 ; n < N ; ++n ){
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
      }
    }
  }
  return(log(denom_res));
}

// [[Rcpp::export]]
NumericVector create_denom_table_bytable_rcpp(NumericVector x, NumericVector mu, NumericVector sigma, int D, int T,  int N, int M){
  if(D * T * N * M != x.length()) stop("D * T * N * M != x.length()");
  NumericVector denom_res(T * N * M);
  denom_res.attr("dim") = IntegerVector::create(T, N, M);
  arma::mat sigma_mat(sigma.length(), sigma.length(), arma::fill::zeros);
  for(std::size_t i = 0 ; i < sigma.length() ; ++i ){
    sigma_mat(i, i) = sigma(i); 
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
  // or(std::size_t i = 0 ; i < T * N * M; ++i ){
  for(std::size_t m = 0 ; m < M ; ++m ){
    for(std::size_t n = 0 ; n < N ; ++n ){
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
      }
    }
  }
  return(log(denom_res));
}

// [[Rcpp::export]]
NumericVector create_weights_table_rcpp(NumericVector loglik_table, NumericVector denom_table, int T,  int N, int M){
  NumericVector weight_table(T * N * M);
  weight_table.attr("dim") = denom_table.attr("dim");
  weight_table = exp(loglik_table - denom_table);
  // weight_table = weight_table/sum(weight_table);
  return(weight_table);
}

// [[Rcpp::export]]
NumericVector gen_mu_chain_apis_rcpp(fp_logposterior_t lp, NumericVector mu, NumericVector sigma, int T, int M){
  int D = mu.length();
  Rcout << "D: "  << D << std::endl ;
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
    xs_curr = gen_xs_rcpp(mu_curr, sigma, D, 1, 1, M);
    loglik_curr = create_loglik_table_rcpp(lp, xs_curr, D, 1, 1, M);
    denom_curr = create_denom_table_bybox_rcpp(xs_curr, mu_curr, sigma, D, 1, 1, M);
    weight_curr = exp(loglik_curr - denom_curr);
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
NumericVector gen_mus_chain_apis_rcpp(fp_logposterior_t lp, NumericMatrix mu, NumericVector sigma, int T, int N, int M){
  int D = mu.ncol();
  NumericVector mu_chains(D * T * N);
  mu_chains.attr("dim") = IntegerVector::create(D, T, N);
  NumericVector mu_curr(D * T);
  IntegerVector idx_res;
  IntegerVector idx_curr;
  for(std::size_t n = 0 ; n < N ; ++n ){
    Rcout << "Je bug !" << std::endl ;
    mu_curr = gen_mu_chain_apis_rcpp(lp, mu(n, _), sigma, T, M);
    for(std::size_t t = 0 ; t < T ; ++t ){
      idx_res = (Rcpp::seq_len(D) - 1) + t * D + n * T * D ;
      Rcout << "idx_res.length() : "  << idx_res.length() << std::endl ;
      idx_curr = (Rcpp::seq_len(D) - 1) + t * D;
      Rcout << "idx_curr.length() : "  << idx_curr.length() << std::endl ;
      mu_chains[idx_res] = mu_curr[idx_curr];
    }
  }
  return(mu_chains);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
logposterior(c(1, 2))
logposterior(c(0, 16))
# fptr1 <- make_GGfuncpp()
# testfun(fptr1, c(1, 5))
lp <- make_logposteterior_cpp()
system.time(mh_rcpp <- simple_MH_rcpp(lp = lp, mu = c(1, 1), sigma = c(3, 3), T = 500))
apply(mh_rcpp, 2, mean)
matplot(mh_rcpp)
T <- 500
N <- 2
sigma <- c(10, 10)
mcmc_chain <- gen_chain_mcmc_rcpp(lp = lp, mu = c(1, 1), sigma = c(3, 3), T = T)
mus_chain <- gen_mus_pmcmc_rcpp(lp = lp, mu = matrix(c(1, 1), ncol = 2, nrow = N), sigma = c(3, 3), T = T, N = N)
mu_chain <- gen_mu_chain_apis_rcpp(lp = lp, mu = c(1, 1), sigma = c(3, 3), T = 500, M = 5)
mus_chain <- gen_mus_chain_apis_rcpp(lp = lp, mu = matrix(c(1, 1), ncol = 2, nrow = N), sigma = c(3, 3), T = 500, N = N,  M = 5)
apply(mus_chain, 1, mean)
xs_chain <- gen_xs_rcpp(mu = mus_chain, sigma = sigma, D = 2, T = T, N = N, M = 2)
loglik_table <- create_loglik_table_rcpp(lp = lp, x = xs_chain, D = 2, T = T, N = N, M = 2)
denom_table_bybox <- create_denom_table_bybox_rcpp(x = xs_chain, mu = mus_chain, sigma = sigma, D = 2, T = T, N = N, M = 2)
denom_table_byrow <- create_denom_table_byrow_rcpp(x = xs_chain, mu = mus_chain, sigma = sigma, D = 2, T = T, N = N, M = 2)
# denom_table_bytable <- create_denom_table_bytable_rcpp(x = xs_chain, mu = mus_chain, sigma = sigma, D = 2, T = T, N = N, M = 2)
weight_table <- create_weights_table_rcpp(loglik_table = loglik_table, denom_table = denom_table_byrow, T = T, N = N, M = 2)
apply(xs_chain, 1, weighted.mean, w = weight_table / sum(weight_table))
library(mvtnorm)
source(file = "R/simple_MH.R")
source(file = "R/luca_posterior.R")
system.time(mh <- simple_MH(lposterior_2, d = 2, T = 5000, mu = c(3, 3), sigma = c(3, 3)))
apply(mh, 2, mean)
matplot(mh)
*/
