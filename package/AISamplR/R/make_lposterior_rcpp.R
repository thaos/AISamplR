#' Creating pointer to a C++ function with Rcpp.
#'
#' \code{make_lposterior_rcpp} helps creating a pointer to a C++ function.
#' 
#' The function helps sourcing the code of a C++ function and
#' returns a pointer to this function that is useable by R. 
#' 
#' To so, the function follows the method describe 
#' \link[=http://gallery.rcpp.org/articles/passing-cpp-function-pointers/]{here}
#' in the Rcpp gallery:
#' http://gallery.rcpp.org/articles/passing-cpp-function-pointers
#' 
#' The function uses \code{RcppArmadillo} and \code{Rcpp}
#' which are imported for the C++ code.
#' Namespace of Rcpp is also imported.
#' 
#' @param body the C++ code of the function that will be sourced by Rcpp. 
#' The function to be exported in R has to be called \code{lposterior}. 
#' \code{lposterior} has to take as argument 
#' an input \code{x} of type Rcpp:: NumericalVector.
#' The function should return a \code{double}.
#' @param ... Additional arguments to be passed to be passed
#' to \code{\link[Rcpp]{sourceCpp}} if necessary.
#' @return A list with the following 2 elements: 
#' \itemize{
##'  \item{fun}{which is the C++ function that imported and useasble in R}
##'  \item{pointer}{an external pointer to the C++ function.
##'   It can be passed as the logposterior argument of 
##'   the functions \code{\link{lais}}, \code{\link{apis}} and \code{\link{pmc}}.}
##' }
#' @examples
#' # Mixture of 2 Gaussian distributions
#' # mu_true = rep(2.5, 8);
#' 
#' lposterior_3 <- function(x){
#'   mu_1 <- c(5, 0)
#'   sigma_1 <- matrix(c(2, rep(0.6, 2), 1), ncol = 2, nrow = 2)
#'   mu_2 <- c(0, 16)
#'   sigma_2 <- matrix(c(3, rep(0, 2), 3), ncol = 2, nrow = 2)
#'   f_1  <- 1/2 * mvtnorm::dmvnorm(x, mean = mu_1 , sigma = sigma_1) 
#'   f_2  <- 1/2 * mvtnorm::dmvnorm(x, mean = mu_2 , sigma = sigma_2) 
#'   f <- f_1 + f_2 
#'   log(f)
#' }
#' 
#' D <- 2
#' T <- 100
#' N <- 50
#' M <- 4
#' 
#' # lais with an R function
#' system.time({
#'   lais_lp3_r <- lais(lposterior_3,
#'       mu = matrix(rnorm(D*N, sd = 10), nrow = D, ncol = N),
#'       sig2_adapt = rep(100, D), sig2_samp = rep(1, D),
#'       compute_logdenom = compute_logdenom_byrow,
#'       N = N, T = T, M = M)
#' })
#' with(lais_lp3_r, compute_expectation(x, weight)) # theorical value: ~ [2.5, 8]
#'
#' #' # With C++ code
#' #' # Create a C++ function and a pointer to the function 
#' body_lp3 <- '
#'   const double log2pi = std::log(2.0 * M_PI);
#'     arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
#'     int n = x.n_rows;
#'     arma::mat x_cen;
#'     x_cen.copy_size(x);
#'     for (int i=0; i < n; i++) {
#'       x_cen.row(i) = x.row(i) - center;
#'     }
#'     return sum((x_cen * cov.i()) % x_cen, 1);    
#'   }
#'   arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false) { 
#'     arma::vec distval = Mahalanobis(x,  mean, sigma);
#'     double logdet = sum(arma::log(arma::eig_sym(sigma)));
#'     arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
#'     
#'     if (log) { 
#'       return(logretval);
#'     } else { 
#'       return(exp(logretval));
#'     }
#'   }
#'   // [[Rcpp::export]]
#'   double lposterior(NumericVector x){
#'     int d = 2;
#'     int d2 = pow(d, 2);
#'     double muarr_1[2] = {5, 0};
#'     double muarr_2[2] = {0, 16};
#'     arma::rowvec mu1(&muarr_1[0], 2);
#'     arma::rowvec mu2(&muarr_2[0], 2);
#'     double sigarr_1[4] = {2, 0.6, 0.6, 1};
#'     double sigarr_2[4] = {3, 0, 0, 3};
#'     arma::mat sigma1(&sigarr_1[0], d, d);
#'     arma::mat sigma2(&sigarr_2[0], d, d);
#'     arma::mat xmat(1, x.length());
#'     for(std::size_t i = 0 ; i < x.length() ; ++i ){
#'       xmat(0, i) = x(i); 
#'     }
#'     double f_1 = sum(dmvnorm_arma(xmat, mu1, sigma1, false));
#'     double f_2 = sum(dmvnorm_arma(xmat, mu2, sigma2, false));
#'     double f = f_1 + f_2;
#'     double ans = log(f) - log(2);
#'     return(ans);
#'   }
#' '  
#' lp3 <- make_lposterior_rcpp(body = body_lp3)
#' 
#' # lais with a pointer to a C++ function
#' system.time({
#'   lais_lp3_rcpp <-  lais(lp3$pointer,
#'     mu = matrix(rnorm(D * N, sd = 3), nrow = D, ncol = N),
#'     sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
#'     compute_logdenom = compute_logdenom_byrow,
#'     N = N, T = T, M = M)
#'  })
#' with(lais_lp3_rcpp, compute_expectation(x, weight)) # theorical value: ~ [2.5, 8]
#' 
#' @export
make_lposterior_rcpp <- function(body, ...){
  head_cpp <- '
    // [[Rcpp::depends(RcppArmadillo)]]
    #include <RcppArmadillo.h>
    using namespace Rcpp;
    typedef double (*fp_logposterior)(NumericVector x);
    typedef XPtr<fp_logposterior> fp_logposterior_t;
  '
  tail_cpp <- '
    // [[Rcpp::export]]
    fp_logposterior_t make_lposterior_cpp(){
      return fp_logposterior_t(new fp_logposterior(lposterior));
    }
  '
  code <- paste(head_cpp, body, tail_cpp, sep = "\n")
  print(code)
  env_cpp <- function(...){
    Rcpp::sourceCpp(code = code, env = environment(), ...)
    list(fun = lposterior, pointer = make_lposterior_cpp())
  }
  env_cpp()
}

