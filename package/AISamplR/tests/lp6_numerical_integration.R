library(devtools)
devtools::load_all()
library(cubature)
library(SimplicialCubature)
# library(rgl)

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