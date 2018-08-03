compute_expectation <- function(xs_chain, weight_table, f = identity, ...){
  weights <- normalize_weights(weight_table)
  xs <- matrix(xs_chain, nrow = dim(xs_chain)[1])
  fxs <- matrix(apply(xs, 2, f, ...), ncol = ncol(xs))
  apply(fxs, 1, weighted.mean, w = weights)
}

rgl_plot <- function(x, y, z){
  ncol  <- 100
  cuts <- cut(z, seq(min(z) , max(z), length.out = ncol + 1))
  col <- rainbow(ncol)[cuts]
  rgl::plot3d(x, y, z)
}

