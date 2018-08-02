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

