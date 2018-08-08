#' AIS: A package for Adaptive Importance Sampling.
#'
#' The AISamplR provides three functions for Adaptive Importance Sampling:
#'   \code{\link{lais}}, \code{\link{apis}} and \code{link{pmc}}.
#'   
#' For a complete description of the three methods,
#' please see the following references:  
#'
#' @docType package
#' @useDynLib AISamplR
#' @name AISamplR
#' 
#' @import Rcpp RcppArmadillo
#' @importFrom grDevices rainbow
#' @importFrom stats dexp dnorm median rmultinom rnorm runif weighted.mean
NULL