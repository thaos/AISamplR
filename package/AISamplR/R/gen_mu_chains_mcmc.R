
gen_mu_chains_mcmc <- function(logposterior,  mu, sigma, T = 100, N = 2){
  D <- ncol(mu)
  if(D != length(sigma)) stop("different dimension for sigma and mu");
  if(N != nrow(mu)) stop("different dimension for N and nrow(mu)");
  mus_res <- array(NA, dim = c(D, T, N))
  for(n in seq.int(N)){
    mus_res[,, n] <- t(simple_MH(logposterior = logposterior,
                                 mu = mu[n, ],
                                 sigma = sigma,
                                 T = T)
    )
  }
  return(mus_res)
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

system.time(
  gen_mu_chains_mcmc(
    lp6$pointer,
    mu = matrix(1:4, ncol = 2, nrow = 100),
    sigma = c(1,1),
    T = 1000, N = 100)
)
system.time(
  gen_mu_chains_mcmc(
    lp6$fun,
    mu = matrix(1:4, ncol = 2, nrow = 100),
    sigma = c(1,1),
    T = 1000, N = 100
  )
)
system.time(
  gen_mu_chains_mcmc(
    lposterior_6,
    mu = matrix(1:4, ncol = 2, nrow = 100),
    sigma = c(1,1),
    T = 1000, N = 100
  )
)

