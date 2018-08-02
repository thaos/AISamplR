library(devtools)
devtools::load_all()
library(cubature)
library(SimplicialCubature)
# library(rgl)

N <- 100
T <- 10000
d <- 2

message(" computing LAIS ...")
set.seed(1)
lais_chain <- lais(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10), sigma = rep(sqrt(3), d), lposterior_6, sigma_mh = rep(1, d), compute_denom = compute_denomtable_byrow)
estim_lais <- compute_expectation(lais_chain$x, lais_chain$w)

message(" computing indep PMC ...")
set.seed(1)
pmc_idchain <- indep_pmc(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_6, compute_denom = compute_denomtable_byrow, reuse_weights = FALSE)
estim_idpmc <- compute_expectation(pmc_idchain$x, pmc_idchain$w)

message(" computing dep PMC ...")
set.seed(1)
pmc_dchain <- pmc(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10) , sigma = rep(sqrt(13), d), lposterior_6, compute_denom = compute_denomtable_byrow, reuse_weights = FALSE)
estim_dpmc <- compute_expectation(pmc_dchain$x, pmc_dchain$w)

message(" computing indep APIS ...")
set.seed(1)
apis_idchain <- indep_apis(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_6, compute_denom = compute_denomtable_byrow)
estim_idapis <- compute_expectation(apis_idchain$x, apis_idchain$w)

message(" computing dep APIS ...")
set.seed(1)
apis_dchain <- apis(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_6, compute_denom = compute_denomtable_byrow)
estim_dapis <- compute_expectation(apis_dchain$x, apis_dchain$w)


message(" Summerizing results ...")
res_table <- matrix(c(estim_lais, estim_idpmc, estim_dpmc, estim_idapis, estim_dapis),
                    nrow = d, ncol = 5)
colnames(res_table) <- c("lais", "idpmc", "dpmc", "idapis", "dapis")


lposterior_6 <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  logtarget <-  -1/32 * (4 - 10 * x1 - x2^2)^2 - x1^2/50 - x2^2/50
}

step <- 0.1
x <- seq(-30, 10, by = step)
y <- seq(-15, 15, step)

xmat <- matrix(x, nrow = length(x), ncol = length(y))
ymat <- matrix(y, nrow = length(x), ncol = length(y), byrow = TRUE)
lp_xy <- apply(cbind(c(xmat), c(ymat)), 1, lposterior_6)
# lp_xy <- apply(cbind(c(xmat), c(ymat)), 1, mvtnorm::dmvnorm)
lp_xy <- exp(lp_xy)
lp_xy <- lp_xy / sum(lp_xy)
lp_xy <- matrix(lp_xy, nrow = length(x), ncol = length(y))
muy <- sum(ymat * lp_xy) #* step^2
mux <- sum(xmat * lp_xy) #* step^2
rgl::plot3d(xmat,ymat, lp_xy, col = rainbow(100)[cut(lp_xy, seq(min(lp_xy), max(lp_xy), length.out = 101))])


# OK?
pi_total <- cubature::adaptIntegrate(function(x){
  exp(lposterior_6(x))
  # mvtnorm::dmvnorm(x)
}, lowerLimit = c(min(x), min(y)), upperLimit = c(max(x), max(y)),
doChecking = TRUE)$integral

mux_adaptIntegrate <- cubature::adaptIntegrate(function(x){
  x[1] * exp(lposterior_6(x)) / pi_total
}, lowerLimit = c(min(x), min(y)), upperLimit = c(max(x), max(y)),
doChecking = TRUE,  tol=1e-10)$integral

# OK
pi_total2 <- integrate(function(y){ 
  sapply(y, function(y) {
    integrate(function(x){
      sapply(x, function(x){
        exp(lposterior_6(c(x, y)))
        # mvtnorm::dmvnorm(c(x, y))
      })}, min(x), max(x), subdivisions = length(x))$value
  })}, min(y), max(y), subdivisions = length(y))$value

mux_doubleintegrate <- integrate(function(y){ 
  sapply(y, function(y) {
    integrate(function(x){
      sapply(x, function(x){
        x  * exp(lposterior_6(c(x, y))) / pi_total2
      })}, min(x), max(x), subdivisions = length(x))$value
  })}, min(y), max(y), subdivisions = length(y))$value

# OK
pi_total22 <- integrate(function(x){ 
  sapply(x, function(x) {
    integrate(function(y){
      sapply(y, function(y){
        exp(lposterior_6(c(x, y)))
        # mvtnorm::dmvnorm(c(x, y))
      })}, min(y), max(y), subdivisions = length(y))$value
  })}, min(x), max(x), subdivisions = length(x))$value

mux_doubleintegrate2 <- integrate(function(x){ 
  sapply(x, function(x) {
    integrate(function(y){
      sapply(y, function(y){
        x  * exp(lposterior_6(c(x, y))) / pi_total22
        # mvtnorm::dmvnorm(c(x, y))
      })}, min(y), max(y), subdivisions = length(y))$value
  })}, min(x), max(x), subdivisions = length(x))$value

#theoretical mu ~= c(0, -1.09)
theo <- c(mu1_theo$integral, mu2_theo$integral)
diff_table <- res_table - theo
rmse <- apply(diff_table, 2, function(x) sqrt(mean(x^2)))

library(rgl)
rgl::plot3d(lais_chain$mu[, 1,], lais_chain$mu[, 2,], 0)
rgl::plot3d(lais_chain$x[, 1,,], lais_chain$x[, 2,,], exp(lais_chain$pi))
rgl::plot3d(lais_chain$x[, 1,,], lais_chain$x[, 2,,], lais_chain$w)



body_lp62 <- '
    double b = 10;
    double sig = 3.5; 
    double x1 = x[0];
    double x2 = x[1];
    double  logtarget = pow(4 - b * x1 - pow(x2, 2)), 2) / (2 * pow(4, 2)) -
      pow(x1, 2) / (2 * pow(sig,2)) - pow(x2, 2) / (2 * pow(sig, 2)); 
    return logtarget;
    }
  ' 

n <- 2
S <- SimplicialCubature::CanonicalSimplex( n )
f0 <- function(x){
  # exp(lposterior_6(x))
  mvtnorm::dmvnorm(x)
}
# Wrong should be equal to one
pi_total3 <- SimplicialCubature::adaptIntegrateSimplex( f0, S )

# Wrong should be equal to one
pracma::integral2(function(x, y){
  # lposterior_6(c(x, y))
  mvtnorm::dmvnorm(x)
},
xmin = min(x), xmax = max(x),
ymin = min(y), ymax = max(y),
vectorized = FALSE)