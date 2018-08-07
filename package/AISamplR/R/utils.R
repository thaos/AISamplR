#' Compute expectation of x
#' 
#' \code{compute_expectation} computes the expectation of the samples \code{x}.
#' 
#' The function computes the expectation of a function of 
#' the samples \code{x} produced by
#' the adaptive importance sampling methods:
#' \code{\link{lais}}, \code{\link{apis}} and \code{link{pmc}}. 
#' If the function \code{f} returns a vector, the expectation is computed
#' for every dimension if the function output.
#' 
#' @param x an array of dimension D x T x M x N with samples drawn
#' from one of the adaptive importance sampling method.
#' @param weight an array of dimension T x M x N with the weight corresponding
#' to every sample \code{x}.
#' @param f a function to be applied to each sample \code{x}.
#' 
#' @examples
#' # draw samples from the "banana shaped distribution" defined by the loglikelihood lposterior_6.
#' lposterior_6 <- function(x){
#'   x1 <- x[1]
#'   x2 <- x[2]
#'   logtarget <-  -1/32 * (4 - 10 * x1 - x2^2)^2 - x1^2/50 - x2^2/50
#' }
#' D <- 2
#' T <- 100
#' N <- 50
#' M <- 4
#' 
#' lais_lp6_r <- lais(lposterior_6,
#'    mu = matrix(rnorm(D * N, sd = 3), nrow = D, ncol = N),
#'    sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
#'    compute_logdenom = compute_logdenom_byrow,
#'    N = N, T = T, M = M)
#' with(lais_lp6_r, compute_expectation(x, weight)) # theorical value: ~ [-1.09, 0]
#' with(lais_lp6_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
#' # expectation of x
#' Ex <- with(lais_lp6_r, compute_expectation(x, weight))
#' #expectation of x^2
#' Ex2 <- with(lais_lp6_r, compute_expectation(x, weight, f = function(x) x^2))
#' # variance of x
#' with(lais_lp6_r, compute_expectation(x, weight, f = function(x) (x - Ex)^2))
#' Ex2 - Ex^2
#' @export
compute_expectation <- function(x, weight, f = identity, ...){
  weight <- normalize_weights(weight)
  xs <- matrix(x, nrow = dim(x)[1])
  fxs <- matrix(apply(xs, 2, f, ...), ncol = ncol(xs))
  apply(fxs, 1, weighted.mean, w = weight)
}



#' 3D scatterplot with rgl.
#' 
#' \code{rgl_plot} makes a 3D scatterplot with a rainbow color scale.
#' 
#' The function is a simple wrapper around \code{\link[rgl]{plot3d}}.
#' It produces a 3D scatterplot with a rainbow color scale for the z-axis.
#' returns a pointer to this function useable by . 
#' 
#' @param x a vector with the x-axis coordinates
#' @param y a vector with the y-axis coordinates
#' @param z a vector with the z-axis coordinates
#' @param ... additional arguments for the function \code{\link[rgl]{plot3d}}.
#' @examples
#' # draw samples from the "banana shaped distribution" defined by the loglikelihood lposterior_6.
#' lposterior_6 <- function(x){
#'   x1 <- x[1]
#'   x2 <- x[2]
#'   logtarget <-  -1/32 * (4 - 10 * x1 - x2^2)^2 - x1^2/50 - x2^2/50
#' }
#' D <- 2
#' T <- 100
#' N <- 50
#' M <- 4
#' 
#' lais_lp6_r <- lais(lposterior_6,
#'    mu = matrix(rnorm(D * N, sd = 3), nrow = D, ncol = N),
#'    sig2_adapt = rep(1, D), sig2_samp = rep(1, D),
#'    compute_logdenom = compute_logdenom_byrow,
#'    N = N, T = T, M = M)
#' with(lais_lp6_r, compute_expectation(x, weight)) # theorical value: ~ [-1.09, 0]
#' with(lais_lp6_r, rgl_plot(x[1,,,], x[2,,,], sqrt(weight)))
#' with(lais_lp6_r, compute_expectation(x, weight))
#' 
#' @export
rgl_plot <- function(x, y, z, ...){
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package \"rgl\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  ncol  <- 100
  cuts <- cut(z, seq(min(z) , max(z), length.out = ncol + 1))
  col <- rainbow(ncol)[cuts]
  rgl::plot3d(x, y, z, col = col, ...)
}

