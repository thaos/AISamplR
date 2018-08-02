#include <RcppArmadillo.h>

typedef double (*fp_logposterior)(Rcpp::NumericVector x);
typedef Rcpp::XPtr<fp_logposterior> fp_logposterior_t;
